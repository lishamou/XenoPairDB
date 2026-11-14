# XenoPairDB data and reproducible analyses

This archive accompanies *XenoPairDB: an evidence-graded resource for organizing pig–human molecular incompatibilities and comparing donor-pig engineering designs* by Lisha Mou, Ying Lu, Zijing Wu, Qishan Zhang and Zuhui Pu.

## Contents

- `analysis/`: frozen scientific inputs, cached public annotations, executable analyses and full-precision expected results. Nested directories preserve the dependency paths used by the analyses; their stage names are organizational identifiers, not alternative biological datasets.
- `figures/`: four manuscript figures (PNG/PDF) with exact source tables and provenance.
- `figure_data/`: numerical source tables for supplementary scorecards, catalogue overlap, sequence QC and contact sensitivity.
- `code/`: deterministic source-table extraction and R plotting scripts.
- `DATA_DICTIONARY.md`, `TABLE_SCHEMAS.tsv`, `ANALYSIS_INDEX.tsv`: interpretation and source locations.
- `provenance/`: source/code crosswalk, source identifiers and reference metadata. `SHA256SUMS.txt` verifies all archived files except itself.
- `LICENSE`, `DATA_LICENSE.md`, `THIRD_PARTY_NOTICES.md`: separate code, original data and upstream reuse conditions.

## Reproduce

Tested scientific environment: Python 3.9.6 with the pinned packages in `requirements.txt`. Create an isolated environment, install these dependencies, then run from the extracted archive:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python reproduce.py --output /absolute/path/to/new_output
```

The output path must not already exist. The replay copies scientific inputs to that path, clears result directories except the explicitly curated sequence-selection input, reruns the computational stages, and compares generated tables with the frozen expected tables. Cached public inputs make scientific replay offline; installation requires package access. No score weights, evidence labels, control definitions or ranking rules are fitted during replay. New runtime/provenance JSON records contain current timestamps and paths and are not expected to be byte-identical to historical records.

The supplied verification report records an executed replay, not a promise of identical behavior under all future software environments. The scientific stages include perturbation, bootstrap and permutation analyses and can take several minutes. Raw mass spectra were not reprocessed, and this archive contains no new wet-laboratory experiment. Historic curation and supplied feature tables are analysis inputs, not outputs claimed to have been recreated from unavailable screening records.

The two display extractors are rerun by `reproduce.py`. To generate graphics from the extracted tables after replay (base R on macOS; Quartz PNG/TIFF support required):

```sh
Rscript code/plot_figures.R /absolute/path/to/new_output
```

The published image exports are supplied from the frozen figure run (R 4.6.2, macOS). The current replay regenerates the four numeric figure source tables; the graphical R export was not rerun in the present environment. Pixel equality is not claimed across graphics devices or font installations. Every figure uses the same frozen numerical data.

## Analysis scope

The catalogue has 75 axes: 59 direct-evidence, 4 disputed and 12 compatibility-control annotations. Sequence-feature comparisons cover 61 axes; stored ISI is numeric for 63. Missing measurements remain missing. Donor masks encode direct genetic targeting, not demonstrated functional rescue or clinical superiority. Score sensitivity and reference-set enrichment are exploratory. Public Liu proteomics and Anand functional data address distinct endpoints and do not independently validate ISI or clinical donor selection. See the article and Supplementary Information for statistical populations, uncertainty and limitations.

## Reuse and citation

Please cite this Zenodo record using the DOI assigned when it is published, and cite the original sources listed in `THIRD_PARTY_NOTICES.md`. Original project code is MIT; original curation and derived materials are CC BY 4.0. Third-party data retain their source licenses, including CC BY-NC 4.0 for the Liu supplements. These third-party files are not relicensed under the project data license. 
