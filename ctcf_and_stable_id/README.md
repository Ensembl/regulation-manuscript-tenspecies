# Finalise the regulatory annotation

These scripts are run after the initial regulatory annotation build. They include code for parsing and adding CTCF features to the base annotation, assigning stable ids and generating the final GFF3 file.

## Prerequisites

The prerequisites are:

- Input data files, which can be downloaded from Zenodo ([https://doi.org/10.5281/zenodo.17496015](https://doi.org/10.5281/zenodo.17496015)):
  - GFF3 files of promoters, enhancers, and unclassified open chromatin regions (base features) from the initial regulatory build i.e., outputs from the regulatory annotation script
  - Some species have CTCF motifs from the Ensembl motif mapping pipeline.
- The `ensreg` Rust package needs to be installed or compiled using Rust. With Rust installed, this can be done with `cargo install ensreg@0.1.0`.
- For Ensembl 114, the scripts were run using:
  - Bash on Ubuntu 24.04.3 LTS
  - bedtools v2.31.1

## Processing steps

### 1. Prepare CTCF features

Some species have CTCF motifs from the Ensembl motif mapping pipeline.

`bedtools` requires `chrom.sizes` files, which are found in the `supporting` directory.

To prepare the CTCF features, run `generate_ctcf_candidate.sh`. For example, for mouse:

```bash
./generate_ctcf_candidate.sh Mus_musculus.GRCm39 ../data/ctcf_and_stable_id/ctcf_motifs/Mus_musculus.GRCm39.CTCF_validated_motifs.bed Mus_musculus_CTCF_candidates.bed
```

### 2. Generate features with stable ids

This step parses the base annotation files and produces the final GFF3 output. With the `-c` option, it includes CTCF features. With `-r`, it replaces the existing ids with stable ids. If mapping from old to new IDs is needed, run the tool prefixed by `RUST_LOG=info`.

For example, for mouse:

```bash
RUST_LOG=info ensreg regbuild mouse -r -c Mus_musculus_CTCF_candidates.bed ../data/ctcf_and_stable_id/base_features/Mus_musculus.GRCm39.regulatory_features.20240615.gff3 > annotation.gff3 2> id_maps.txt
```

### 3. Sort the annotation

This can be done with `bedtools sort -header`. For example:

```bash
bedtools sort -header -i annotation.gff3 > annotation_sorted.gff3
```
