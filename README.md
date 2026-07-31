# MOCAP: MOtif Characterization & Annotation Pipeline
MOCAP is a comprehensive and configurable bioinformatics pipeline designed to analyze DNA motifs discovered from deep learning models. It automates the entire workflow, from extracting raw motifs from modisco HDF5 files to generating a final, interactive report for key genes.

The pipeline is built in Python and designed to be modular, memory-efficient, and easily configurable through a central config.yaml file. It replaces a series of R scripts with a single, streamlined application.

## Features
**Configurable Workflow:** Control the entire pipeline from a single, human-readable .yaml file.

**Modular Steps:** Each stage of the analysis is a separate, well-defined step.

**Multiple Motif Types:** Extract and process Position Frequency Matrices (PFMs), Position Weight Matrices (PWMs), or Contribution Weight Matrices (CWMs).

**Motif Cleaning:** Optional trimming of noisy, low-information content positions from motif ends.

**Comprehensive Comparison:** Compare your motifs against multiple, curated public databases (JASPAR, Arabidopsis DAP, RNA-binding proteins, miRNAs) with separate, organized outputs.

**Positional Analysis:** Calculate detailed statistics on the genomic location of motifs relative to transcription start and end sites (TSS/TTS).

**Genomic Annotation:** Map motif occurrences to a reference genome and annotate them based on a GFF/GTF file.

**Advanced Visualization:** Generate high-quality sequence logos for all motifs and smoothed positional density plots. Produce interactive, JBrowse-like HTML reports for specific genes, complete with a detailed summary table.

**Memory Efficient:** Designed to handle large occurrence files by processing data in manageable chunks.

**Caching System:** Automatically skips steps that have already been completed to save time on re-runs.

##  Project Structure
The project is organized into a clean, standard Python package structure.

moca/
```
├── configs/              # All .yaml configuration files and requirements.txt
├── data/                 # Input data (HDF5 models, reference genomes, etc.)
├── results/              # All output files, organized by run name
├── scripts/              # Helper shell scripts (e.g., for blamm)
├── run_pipeline.py       # The main script to execute the entire pipeline
├── setup.py              # Makes the project installable
└── moca/                 # The core Python source code package
    ├── __init__.py
    └── steps/            # Each pipeline step is a separate script
        ├── __init__.py
        ├── step1_nomenclature.py
        └── ... (all other step scripts)
```

## Setup and Installation
The pipeline uses a Conda environment to manage its dependencies.

**Create the Conda Environment:**

`conda create --name moca-env python=3.9`

**Activate the Environment:**

`conda activate moca-env`

**Install Dependencies:** 
The configs/ directory should contain a requirements.txt file. Install all necessary libraries with a single command:

`pip install -r configs/requirements.txt`

**Install MOCA:** 
To make the pipeline's code accessible to Python, install it in "editable" mode from the main project directory:

`pip install -e .`

### External requirements
MOCAP was optimized and tested for the modisco output of DeepCRE (https://github.com/NAMlab/deepCRE_reimplemented/).  
MOCAP used BLAMM (https://github.com/biointec/blamm) for motif search in fasta-sequences like genomic data. Please install BLAMM and set the PATH in the config.yaml file and respective scripts-file. 
MOCAP optimizes the search of genomic sequence via partitioning. MOCAP was optimized and tested to work with samtools (https://www.htslib.org/). Please install samtools, or others (e.g. agat; https://agat.readthedocs.io/en/latest/) and set the PATH in the config.yaml file and respective scripts-file.
MOCAP uses convenient publicly available input data in addition to the modisco feature extraction. Please set up the following directory with recommended files: 
*motif_database* from https://meme-suite.org/meme/doc/download.html
*reference* this directory shall contain species specific genome (.fa, .fas) and annotation files (.gff, .gff3, .gtf)
*deepCRE_results* this directory shall contain the results of https://github.com/NAMlab/deepCRE_reimplemented, modisco results 

`cd data`
`mkdir motif_database reference deepCRE_results`
`cd ..`

## Workflow
Running the pipeline is a simple, two-step process: configure and execute.


### 1. Configure Your Analysis
All aspects of a pipeline run are controlled by a .yaml file inside the configs/ directory. You can create a new file for each analysis (e.g., configs/arabidopsis_run.yaml).

**Key Configuration Parameters:**

run_name: A unique name for your analysis. All results will be saved in a folder with this name inside results/.

input_hdf5: The path to the modisco HDF5 file you want to analyze.

species_tag & model_tag: Short identifiers for your species and model, used for naming output files.

use_trimmed_motifs: A global switch (true or false) that tells all downstream steps whether to use the raw or trimmed motif files.

Step Sections: Each step of the pipeline has its own section (e.g., nomenclature, comparison, annotation). To run a step, set run: true within its section.

### 2. Execute the Pipeline
Once your configuration file is ready, run the entire pipeline from your main moca/ directory with a single command:

`python run_pipeline.py --config configs/your_config_file.yaml`

The pipeline will execute each step marked with run: true, printing its progress to the terminal.

### Pipeline Steps Explained
Nomenclature: Extracts motifs from the HDF5 file. Can be configured to extract PFMs (for sequence patterns) or CWMs (for contribution scores). Optionally trims noisy ends from motifs.

Visualization: Generates high-quality PNG sequence logos for each motif.

Comparison: Compares your motifs against a list of local public databases and generates separate reports for each.

Importance: Calculates metadata for each motif, such as the sum of its contribution scores.

Saliency: (Optional) Visualizes saliency maps for specific genes.

Clustering: Performs hierarchical clustering of your motifs and generates a dendrogram.

Ranging: Analyzes the positional distribution of motifs and calculates statistics relative to the TSS and TTS. Places both halves of the model input on one axis and records the frame in `<spec><model>-ranging_geometry.json`.

Projection: Maps your motifs to a reference genome using the external blamm tool (https://github.com/biointec/blamm).

Annotation: Merges the motif occurrences with a GFF/GTF file to annotate which motifs fall within which gene regions, and applies the per-EPM positional filters. See "Positional filtering".

Browser Viz: Generates a comprehensive, interactive HTML report for specific genes, combining a JBrowse-like gene model with a detailed summary table.

Performance: (Optional) Evaluates the predictive performance of your motifs against gene expression data.

## Positional filtering

The positional filter is **per EPM and per region** — each EPM has its own band, learned from where
deepCRE recognises it. There is no shared window: on a 48-EPM Arabidopsis model the q10–q90 bands range
from 150 to 831 bp wide, and no single position in the window is accepted by all EPMs.

Ranging places both halves of the model input on one axis: `1` = outer edge of the flank,
`upstream_len + 1` = the gene border (TSS or TTS), up to `upstream_len + body_len` = `body_len` into the
transcript. Annotation converts each BLAMM hit into that same frame, so thresholds and hits are always
directly comparable. Ranging writes the frame to `<spec><model>-ranging_geometry.json`, and annotation
reads it rather than assuming it.

Window geometry is configurable for non-vanilla deepCRE setups and for phytoExpr. Defaults are deepCRE
vanilla — 1000 bp upstream of the TSS + 500 bp into the transcript, a 20 bp N spacer, then 500 bp from the
transcript 3' end + 1000 bp downstream of the TTS, totalling 3020 bp. Existing configs that omit the
`window:` and `annotation.filters:` blocks keep working on these defaults.

```yaml
window:
  upstream_len: 1000
  body_len: 500
  spacer: 20

annotation:
  dedupe: true              # collapse the 2x F/R redundancy, keyed on (gene, position, epm)
  deep_intragenic: drop     # drop | pass  — hits deeper than body_len from either border
  filters:
    band: q1q9              # none | minmax | q1q9
    weight_region:
      enabled: true
      min_region_fraction: 0.15
    min_region_seqlets: 20
  strict_merge:
    min_match_rate: 0.5
```

**`band`** is the per-EPM plausibility test. Note it is *shape-adaptive* rather than a dispersion filter:
the inner 80% of a broad seqlet distribution is a broad band, so retention scales with band width
(Spearman 0.72 — the narrowest quartile of bands retains ~13% of in-scope hits, the widest ~45%). Broad
bands are genuine signal, reflecting EPMs deepCRE recognises over a wide region, so they are not penalised.

**`weight_region`** is the abundance gate, and the one that removes poorly supported EPMs. Across five runs
it lifted the minimum retained band width from 0 bp to 120–143 bp while leaving the median and the broad
tail untouched, and raised median seqlet support by 15–40%, at a cost of 2.5–4.3 percentage points of
occurrences. A band built from a single seqlet has `q10 == q90`, i.e. a width of one base pair; these are
what it removes.

**`min_region_seqlets`** is an absolute floor, since `weight_region` is relative and misses EPMs that are
rare in both regions.

**`strict_merge.min_match_rate`** aborts the run if too few in-range occurrences match a gene key from
`reference_gff`, guarding against silently annotating a fraction of the data when the GFF build does not
match the one the projection was scanned against.

`iqr` and `sd` remain in the ranging output as descriptive columns but are deliberately not filters. `iqr`
*is* `q90 - q10`, so filtering on it removes the broadest bands — the opposite of the intent — and leaves
the degenerate zero-width bands untouched.

Annotation writes `annotation_report.json` with the geometry, filter settings, GFF match rate, the full
filter cascade and every dropped band, so a run stays auditable afterwards.

## Corrections in the positional filter

The annotation step previously computed `dist_transc_border` as a distance measured **outward from the gene
border**, domain `[0, flank_size]`, and compared it directly against `q10`/`q90` — which ranging expresses
in a frame where the gene border sits at ~1000 and each EPM's preferred range typically runs 1000–1490,
i.e. **inside** the transcript. The two frames are mirrored and offset.

Measured across six runs (three model tags, five target genomes, self- and cross-species projection):

- `q10` exceeded the flank width for 20–36 of every ~50 EPMs per region, making their bands unreachable
- **45–60 of 58–111 EPM×region pairs retained exactly zero hits** — annihilated, not thinned
- 88–97% of hits, everything classified `intragenic`, bypassed the filter entirely
- flank retention was 1.3–9.9%, against 18.0–49.6% once the frame is corrected

Annotation now places each hit into the ranging frame rather than converting thresholds, which is also
robust to chromosome-edge truncation of the extracted region. For the record the fault was **not** a
missing `-1520` TTS offset: that offset belongs to an older ranging convention
(`moca_blue/mo_ran/..._TSS-TTS.1.4.R`) and is correctly commented out in the modern reference scripts.
Strand handling was correct throughout — gene strand and BLAMM strand are symmetric both in scope and
after filtering in every run tested.

Ranging fixes:

- the TTS half was flipped with `sequence_length - start`, placing the TTS border at 1000 while the TSS
  border sat at 1001; now `(sequence_length + 1) - start`, so both halves agree
- the TTS range bound was `[1520, 3000]`, admitting one base of the N spacer and discarding the 20 most
  distal downstream bases while the TSS half used its full 1500; now `[1521, 3020]`
- `mode` is binned after the flip, so it shares the frame of every other column
- `cv` is NaN-guarded (`sd` is undefined for single-seqlet patterns)

Restored from the moca_blue reference (`mo_proj/mo_feat-filter.v3.x`), all config-driven: nearest-border
assignment so no hit is assigned to both borders in a short gene; the in-scope pre-filter, now applied per
chunk *before* the GFF merge, which also fixes the previous concat-then-filter memory profile;
deduplication; and `weight_region`.

The legacy `word_size` gate was deliberately **not** restored. BLAMM emits exactly one span per EPM variant
(verified: 0 of 100 and 0 of 112 EPMs across two runs show more than one), so it never removed partial
matches — its real effect was rejecting border-straddling hits, which nearest-border assignment now
handles directly.

Further fixes: `flank_size` is read from config instead of a hardcoded constant, and the run aborts if it
is smaller than `window.upstream_len` (a config setting 1500 previously produced a silent empty result,
because the merge key could never match); the gene key reproduces the projection step's
`max(1, start - flank)` clamp, so genes near a chromosome start no longer lose hits silently; dedupe
orders the forward CWM first so the surviving `epm_instance` and `strand` are deterministic; BED output is
0-based half-open (`chromStart` was previously the 1-based coordinate); and `epm` is the canonical
identifier column, with `epm_instance` carrying the F/R suffix and `cwm_strand` recording which CWM
orientation matched.
