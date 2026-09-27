# Single-cell multi-omics reveals cell-type-resolved regulatory programs and genetic prioritization in lupus nephritis

This repository contains computational scripts used for the integrative analysis of single-cell transcriptomics, single-cell chromatin accessibility, genetic variation, and regulatory genomics in lupus nephritis (LN).

The workflow integrates:

- single-cell RNA sequencing (scRNA-seq)
- single-nucleus ATAC sequencing (snATAC-seq)
- peak-to-gene regulatory linkage analysis
- transcription factor (TF) regulatory network inference
- lupus nephritis blood eQTL mapping
- GWAS fine-mapping and SNP heritability analysis
- SNP-to-cis-regulatory-element (CRE) mapping
- sequence-based regulatory prediction using gkm-SVM

The overall analytical framework is designed to identify cell-type-specific regulatory programs and prioritize candidate genes associated with systemic lupus erythematosus (SLE) and kidney function traits.

---

# Overview of computational workflow
             Lupus nephritis multi-omics integration

                        scRNA-seq
                           |
                           |
          Immune cell states and transcriptional programs
                           |
                           |
                        snATAC-seq
                           |
                           |
      Chromatin accessibility, DARs, TF activity and motifs
                           |
                           |
                    Peak-to-gene linkage
                           |
                           |
             TF–CRE–gene regulatory networks
                           |
          -------------------------------------
          |                                   |
        eQTL                             GWAS
          |                                   |
          -------------------------------------
                           |
                           |
         Colocalization and genetic prioritization
                           |
                           |
          Candidate regulatory genes and pathways


