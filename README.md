# Lineage-aware analysis of AprX-associated proteolytic spoilage

Code accompanying the manuscript:

**Lineage structure limits sequence-based prediction of AprX-associated proteolytic spoilage in Pseudomonas**

## Contents

The python scripts in the repository contains the staged analysis workflow used for:

- phenotype extraction and integration
- AprX/serralysin sequence screening
- SVMSY-region mapping
- core-genome phylogenetic analysis
- lineage-aware phenotype analysis
- external RefSeq type-material validation
- serralysin reference-proximity analysis
- aprX-lipA2 locus reconstruction
- operon architecture analysis
- independent zinc-anchored validation of the serralysin Met-turn region
- validated SVMSY statistical reanalysis
- generation of publication summaries and figures

## Data

This is a code-only repository.

Genome sequences, intermediate datasets, result tables, figures and manuscript files are not distributed here.

Public genome accessions and source datasets are described in the manuscript and supplementary material.

## Software

The workflow uses Python together with external bioinformatics tools including:

- BLAST+
- MAFFT
- IQ-TREE 2
- HMMER

Some Python scripts also use packages including pandas, NumPy, SciPy, matplotlib and Biopython.

## Reproducibility

Scripts are numbered approximately in workflow order. Intermediate files are generated during the workflow and are intentionally excluded from this repository.
