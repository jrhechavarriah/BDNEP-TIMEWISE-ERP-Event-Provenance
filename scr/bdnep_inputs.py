#!/usr/bin/env python3
"""Verified replay artifact locator for BDNEP repository or original P14 directory.

Supported layouts:
  <repo>/SHA256SUMS.txt and <repo>/results/<name>
  <P14 replay>/SHA256SUMS.txt and <P14 replay>/<name>
No unverified fallback; no EEG processing; no file writes.
"""
from __future__ import annotations
import hashlib
import re
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verified_artifact(source: Path, name: str) -> Path:
    source = Path(source).resolve()
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', name) or name in ('.', '..'):
        raise RuntimeError('Unsafe artifact name: ' + name)
    manifest = source / 'SHA256SUMS.txt'
    if not source.is_dir() or not manifest.is_file() or manifest.is_symlink():
        raise RuntimeError('Missing verified source root/SHA256SUMS.txt: ' + str(source))
    repo_layout = (source / 'results').is_dir() and (source / 'reference').is_dir()
    rel = ('results/' if repo_layout else '') + name
    path = source.joinpath(*rel.split('/'))
    if not path.is_file() or path.is_symlink():
        raise RuntimeError('Required replay artifact missing: ' + rel)
    matches = []
    seen = set()
    for idx, raw in enumerate(manifest.read_text(encoding='utf-8-sig').splitlines(), 1):
        if not raw.strip():
            continue
        parts = raw.split(maxsplit=1)
        if len(parts) != 2 or not re.fullmatch(r'[0-9a-fA-F]{64}', parts[0]):
            raise RuntimeError(f'Malformed SHA256 manifest line {idx}')
        entry = parts[1].strip().lstrip('*').replace('\\', '/')
        entry_parts = entry.split('/')
        if entry in seen or any(x in ('', '.', '..') for x in entry_parts) or ':' in entry or entry.startswith('/'):
            raise RuntimeError('Unsafe or duplicate SHA256 entry: ' + entry)
        seen.add(entry)
        if entry == rel:
            matches.append(parts[0].lower())
    if len(matches) != 1 or sha256(path) != matches[0]:
        raise RuntimeError('Replay input SHA256 mismatch: ' + rel)
    return path
