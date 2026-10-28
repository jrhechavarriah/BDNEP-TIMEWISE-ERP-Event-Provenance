#!/usr/bin/env python3
"""Acquire/verify SHA256-pinned NEMAR and OSF data OUTSIDE the BDNEP repository.

Modes:
  --manifest <repo>/reference/source_manifest_v1.json --destination <EXTERNAL_DATA_ROOT>
  --manifest ... --nemar-dir <EXTERNAL_NEMAR_DIR> --osf-dir <EXTERNAL_OSF_DIR>
Add --verify-only to check existing data without network or writes.
All 261 URLs, relative paths and target paths are validated before downloading.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from urllib.error import HTTPError, URLError

PINNED_MANIFEST_SHA = 'f75505d2b30f57d4b6ce88a5efcd60c2cf27398e7021dbe2545f34dab4d4827c'
CHUNK = 1024 * 1024


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(CHUNK), b''):
            result.update(chunk)
    return result.hexdigest()


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def inside(path: Path, parent: Path) -> bool:
    path = path.resolve(strict=False)
    parent = parent.resolve(strict=False)
    return path == parent or parent in path.parents


def prepare(manifest_path: Path, destination: Path | None, nemar_dir: Path | None,
            osf_dir: Path | None):
    manifest_path = manifest_path.resolve(strict=True)
    require(digest(manifest_path) == PINNED_MANIFEST_SHA, 'Pinned manifest SHA256 mismatch')
    data = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    require(data.get('file_count') == 261 and len(data.get('files', [])) == 261,
            'Incomplete or modified source manifest')
    # Manifest is located inside the BDNEP/reference directory.
    require(manifest_path.parent.name.lower() == 'reference',
            'Manifest must be inside <BDNEP>/reference to enforce output isolation')
    repo = manifest_path.parent.parent.resolve()
    require((repo / 'scr').is_dir(), 'Not a BDNEP repository layout')
    if destination is not None:
        require(nemar_dir is None and osf_dir is None,
                'Use --destination OR both --nemar-dir/--osf-dir')
        dest = destination.resolve(strict=False)
        roots = {'nemar': dest / 'nemar', 'osf': dest / 'osf'}
    else:
        require(nemar_dir is not None and osf_dir is not None,
                'Supply --destination OR both --nemar-dir and --osf-dir')
        roots = {'nemar': nemar_dir.resolve(strict=False), 'osf': osf_dir.resolve(strict=False)}
    for label, folder in roots.items():
        require(not inside(folder, repo), 'Refusing source data inside repository: ' + label)
        require(not inside(folder, manifest_path.parent), 'Source root overlaps reference dir')
    files = []
    seen = set()
    for entry in data['files']:
        rel = entry['path'].replace('\\', '/')
        parts = rel.split('/')
        require(len(parts) >= 2 and parts[0] in roots and all(x not in ('', '.', '..') for x in parts)
                and ':' not in rel and not rel.startswith('/'), 'Unsafe source path: ' + rel)
        require(rel not in seen, 'Duplicate source path: ' + rel)
        seen.add(rel)
        url = entry['url']
        parsed = urllib.parse.urlsplit(url)
        require(parsed.scheme == 'https' and parsed.username is None and parsed.password is None,
                'Unapproved scheme or credentials: ' + rel)
        if parts[0] == 'nemar':
            require(url.startswith('https://data.nemar.org/nm000132/v1.1.1/'),
                    'NEMAR URL outside pinned dataset/version')
        else:
            require(re.fullmatch(r'https://osf\.io/download/[a-zA-Z0-9]+/?', url) is not None,
                    'Noncanonical OSF URL: ' + rel)
        expected = str(entry['sha256']).lower()
        require(re.fullmatch(r'[a-f0-9]{64}', expected) is not None,
                'Missing source hash: ' + rel)
        path = roots[parts[0]].joinpath(*parts[1:])
        require(not inside(path, repo), 'Destination resolves into repository: ' + rel)
        require(inside(path, roots[parts[0]]), 'Destination escapes source root: ' + rel)
        files.append((rel, url, path, expected))
    return files


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', required=True, type=Path)
    dst = p.add_mutually_exclusive_group()
    dst.add_argument('--destination', type=Path)
    dst.add_argument('--nemar-dir', type=Path)
    p.add_argument('--osf-dir', type=Path)
    p.add_argument('--verify-only', action='store_true')
    args = p.parse_args()
    files = prepare(args.manifest, args.destination, args.nemar_dir, args.osf_dir)
    count = 0
    created = 0
    for idx, (rel, url, path, expected) in enumerate(files, 1):
        # Recheck resolution to reject newly introduced symlinks before filesystem access.
        if path.is_file():
            require(digest(path) == expected, 'Existing source SHA256 differs: ' + rel)
            print(f'[{idx}/261] VERIFIED EXISTING {rel}', flush=True)
        elif args.verify_only:
            raise RuntimeError('Missing source in --verify-only mode: ' + rel)
        else:
            require(not path.is_symlink(), 'Unsafe output symlink: ' + rel)
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(path.name + '.bdnep.part')
            require(not tmp.is_symlink(), 'Unsafe temporary file symlink: ' + rel)
            last = None
            for attempt in range(3):
                try:
                    req = urllib.request.Request(url, headers={'User-Agent': 'BDNEP-TIMEWISE-Reproducibility/1.0'})
                    with urllib.request.urlopen(req, timeout=120) as response, tmp.open('wb') as out:
                        for chunk in iter(lambda: response.read(CHUNK), b''):
                            out.write(chunk)
                    require(digest(tmp) == expected, 'Downloaded SHA256 differs: ' + rel)
                    tmp.replace(path)
                    last = None
                    break
                except (OSError, URLError, HTTPError, RuntimeError) as exc:
                    last = exc
                    if attempt < 2:
                        time.sleep(attempt + 2)
            if last is not None:
                raise RuntimeError(f'Source download failed: {rel}: {last}')
            created += 1
            print(f'[{idx}/261] VERIFIED DOWNLOAD {rel}', flush=True)
        count += 1
    print('SOURCE_FILES_SHA256_VERIFIED:', count, '/261')
    print('NEWLY_DOWNLOADED:', created)


if __name__ == '__main__':
    main()
