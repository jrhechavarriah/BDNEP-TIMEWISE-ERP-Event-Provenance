# BDNEP-TIMEWISE — Reproducible ERP event-marker provenance analysis

Version **1.0.0**. Software author: **Jesus Rafael Hechavarria-Hernandez** (ORCID 0000-0002-9013-8665).
Software DOI (reserved until Zenodo publication): **10.5281/zenodo.23567028**.

## Overview
Independent replay of ERP CORE N170 data, combining the NEMAR nm000132 v1.1.1 continuous EEG with OSF original and author-corrected event labels. The documented digital-marker difference is +26 samples for stimuli at 1024 Hz; physical screen onset was **not** calibrated. The primary contrast is PO8 face-minus-car amplitude in the prespecified 110–150 ms window, Arm A. Arm B and PO7 results are descriptive. No artifact rejection, ICA or bad-channel interpolation.

## Data and rights
EEG and OSF event-source bytes are **NOT** redistributed in this repository. Download them from their originating NEMAR and ERP CORE/OSF sources, subject to source-specific licenses and required attribution. Original BDNEP-TIMEWISE Python software is licensed under GPL-3.0-only as specified in LICENSE. This license does not apply to upstream datasets, event-source files, or independently licensed dependencies. Third-party libraries retain their own licenses.

## Reproduce (run commands from repository root)

1. Install Python dependencies with `python -m pip install -r requirements.txt` (the pinned environment is documented separately).
2. Download the 261 source files **outside** the repository:
   `python scr/00_download_sources.py --manifest reference/source_manifest_v1.json --destination EXTERNAL_DATA_DIRECTORY`
3. Verify files:
   `python scr/00_download_sources.py --manifest reference/source_manifest_v1.json --destination EXTERNAL_DATA_DIRECTORY --verify-only`
4. Independent full cohort replay (requires EEG files):
   `python scr/02_replay_erp.py --reference reference --nemar EXTERNAL_DATA_DIRECTORY/nemar --osf EXTERNAL_DATA_DIRECTORY/osf --output-root EXTERNAL_OUTPUT_DIRECTORY --all`
5. Rebuild tabulated results and main/supplementary figures from **verified** replay results or the supplied `results/` snapshot:
   `python scr/07_generate_tables.py --repo . --out-dir EXTERNAL_OUTPUT_DIRECTORY/tables`
   `python scr/03_generate_figure2.py --repo . --out-dir EXTERNAL_OUTPUT_DIRECTORY/figure2`
   `python scr/04_generate_supplementary_S1.py --repo . --out-dir EXTERNAL_OUTPUT_DIRECTORY/s1`
   `python scr/06_generate_figure3_po8_arma.py --repo . --out-dir EXTERNAL_OUTPUT_DIRECTORY/figure3`
6. Tests:
   `python -B -m unittest discover -s scr -p test_supplementary_generators.py -v`
   `python scr/02_replay_erp.py --reference reference --verify-contract`

Figures S2 (Face/Car behavioral eligibility) and S3 (descriptive unfiltered EEG anomalies) are identified in `supplementary/FIGURE_NUMBERING.md`.

## Reproducibility notes
The public `reference/` directory contains scientific protocol, cohort, and source-provenance records for the reproducible distribution. Administrative release-gate fields from the historical internal workflow are excluded from the public package. The SHA-256 values in `SHA256SUMS.txt` describe this public release's byte-identical files; they do not replace the separately archived historical source and replay records.

The `results/` directory contains the independently reproduced numerical outputs, not the original EEG source recordings. Public scientific reproducibility does not imply permission to redistribute upstream raw data. See `CITATION.cff` for the software citation and cite the original ERP CORE authors/data providers separately.
