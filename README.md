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
---

# Repository structure


.
├── scRNA-seq/
├── snATAC-seq/
├── TF-regulators/
├── Peak2Gene Analysis/
├── eQTL Analysis/
├── Genotype_data_analysis/
├── SNP-CRE/
├── SNP_heritability_Analysis/
├── gkm-SVM/
└── RNA-seq/


---

# 1. scRNA-seq analysis

Directory:


scRNA-seq/


This folder contains scripts for processing and analyzing PBMC single-cell RNA sequencing data.

Main analyses include:

- quality control and visualization
- dataset integration
- clustering and cell-type annotation
- immune state characterization
- differential expression analysis
- pseudobulk differential expression
- cell-cell communication analysis
- CIBERSORTx reference generation

Main scripts:

| Script | Description |
|---|---|
| `1.scRNA-seq-QC-violinPlot.r` | QC visualization |
| `2.Final10hc_11LN_integration.R` | scRNA-seq integration and downstream analysis |
| `CellChat.r` | Cell-cell communication analysis |
| `ExpressionMatrix_Reference_forCIBERSORT.r` | Generate reference matrix for immune fraction estimation |
| `prepare-For-sccoda.r` | Prepare data for compositional analysis |

Pseudobulk analysis:


scRNA-seq/pseudobulk/


Includes:

- cell-type-specific DEG analysis
- edgeR-based pseudobulk testing
- DESeq2 comparison
- sensitivity analysis excluding highly treated samples


---

# 2. snATAC-seq analysis

Directory:


snATAC-seq/


This folder contains scripts for chromatin accessibility analysis using ArchR.

Analyses include:

- quality control
- dimensional reduction
- clustering
- cell-type annotation
- differential accessibility analysis
- peak calling
- TF motif analysis
- regulatory element identification


Main scripts:

| Script | Description |
|-|-|
| `QC.r` | Quality control analysis |
| `ATAC-processing.r` | ArchR preprocessing workflow |
| `Subcluster_analysis.r` | Cell-type-specific snATAC analysis |
| `ArchR_HRG_peakClustering.r` | High regulatory gene analysis |


---

# 3. Peak-to-gene regulatory linkage analysis

Directory:


Peak2Gene Analysis/


and


snATAC-seq/Permutation/


Peak-to-gene linkage analysis was performed using ArchR and extended with permutation-based significance evaluation.

Analyses include:

- ArchR peak-to-gene linkage inference
- permutation-based peak-to-gene significance estimation
- overlap-based sensitivity analysis
- regulatory peak clustering
- high regulatory gene (HRG) identification


Important scripts:

| Script | Description |
|-|-|
| `CallPeak_addPeak2GeneLinks.r` | Peak calling and ArchR peak-to-gene linkage |
| `ArchR_addPeak2GeneLinks_subcluster.r` | Cell-type-specific peak-to-gene linkage |
| `addPermPeak2GeneLinks.R` | Permutation-based peak-to-gene testing |
| `AddPermPeak2gene_getp2GR.r` | Generate permutation-adjusted regulatory links |
| `addPerm_Overlapcutoff_358_SensitivityAnalysis.r` | Sensitivity analysis |
| `ArchR_peakClustering_PermOC5_pseudobulk.r` | Regulatory peak clustering |

---

# 4. Leave-one-donor-out validation of peak-to-gene links

Directory:


snATAC-seq/LODO_P2G_scripts/


This folder contains donor-level robustness analysis.

Purpose:

To evaluate the reproducibility of peak-to-gene regulatory links by removing one LN donor at a time.

Workflow:


Full dataset peak-to-gene links

      |
      |

Leave one donor out

      |
      |

Recalculate regulatory links

      |
      |

Compare overlap and reproducibility



Main scripts:

| Script | Description |
|-|-|
| `01_LODO_setup.R` | Prepare leave-one-donor-out datasets |
| `03_LODO_P2G.R` | Calculate LODO peak-to-gene links |
| `04_compare_full_LODO.R` | Compare full and LODO results |
| `05_HRG_reproducibility.R` | Evaluate HRG robustness |

---

# 5. Transcription factor regulatory analysis

Directory:


TF-regulators/


This folder contains scripts for identifying TF-driven regulatory programs.

Analyses include:

- TF activity estimation
- TF target identification
- TF-associated CRE analysis
- cell-type-specific regulatory programs


Main script:


TF_regulators_identification.r



---

# 6. eQTL analysis

Directory:


eQTL Analysis/


This folder contains scripts for blood eQTL mapping in lupus nephritis samples.

Analyses include:

- genotype-expression preprocessing
- covariate selection
- PEER factor optimization
- MatrixeQTL association testing



Main scripts:

| Script | Description |
|-|-|
| `MatrixEQTL_10peer_nokgp.R` | eQTL mapping using MatrixeQTL |
| `MatrixEQTL_10peer_nokgp_noCF.R` | eQTL analysis without cell fraction adjustment |
| `coloc.r` | Colocalization analysis |
| `find_final_eqtl.R` | Identify final significant eQTL pairs |

---

# 7. Genotype quality control and imputation

Directory:


Genotype_data_analysis/


Contains scripts for genotype preprocessing.

Includes:

- genotype quality control
- genotype imputation


Scripts:


Quanlity_Control.sh
Imputation.sh



---

# 8. SNP-CRE regulatory annotation

Directory:


SNP-CRE/


This folder contains scripts linking fine-mapped variants to candidate cis-regulatory elements.

Main analysis:


Links Fine-mapped SNPs to CRE.r



Workflow:


Fine-mapped SNPs

   |

Chromatin accessible regions

   |

Candidate CREs

   |

Target genes



---

# 9. GWAS fine-mapping and SNP heritability analysis

Directory:


SNP_heritability_Analysis/


Includes:

- GWAS annotation
- SNP heritability enrichment analysis
- LDSC analysis
- fine-mapped variant interpretation


Scripts:

| Script | Description |
|-|-|
| `LDSC.sh` | LD score regression |
| `FinemappedSNP-analysis.r` | Fine-mapped SNP analysis |


---

# 10. gkm-SVM regulatory sequence modeling

Directory:


gkm-SVM/


This folder contains sequence-based prediction analysis for regulatory elements.

Includes:

- SNP sequence extraction
- negative sequence generation
- gkm-SVM training
- prediction performance evaluation


Scripts:

| Script | Description |
|-|-|
| `get_snp_seqs.py` | Extract SNP sequences |
| `get_null_seqs.R` | Generate background sequences |
| `gkmsvm_predict_CV.pl` | Cross-validation prediction |
| `compute_gkmsvm_cv_performance.r` | Model performance evaluation |

---

# Software requirements

Main computational environments:

- R >= 4.1
- Python >= 3.8
- ArchR
- Seurat
- edgeR
- DESeq2
- MatrixEQTL
- coloc
- LDSC
- gkm-SVM


Specific package versions should be considered when reproducing the analyses.


---


