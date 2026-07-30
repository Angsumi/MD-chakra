# Antigravity (AGY) Context & Project History

This file serves as a memory core for AGY to understand the transcriptomic analysis project structure, the analytical decisions made, and the current state of the workspace.

## 1. Project Overview
*   **Project Name:** MD Chakra Transcriptomics
*   **Samples:** 4 individual samples (`NPFM1`, `NPFM2`, `NPFM3`, `NPFM4`).
*   **Crucial Analytical Context:** **There are NO biological replicates and NO experimental groups (e.g., Control vs. Treatment).** 
    *   An earlier iteration of the analysis erroneously grouped these into artificial "Control" and "Treatment" groups, running standard DESeq2 statistical differential expression (DE).
    *   **Action Taken:** All invalid DE tables, MA plots, and Volcano plots resulting from that false assumption were permanently deleted. The analysis was completely rewritten to treat all 4 samples completely independently.

## 2. Preprocessing & Upstream Analysis
The upstream pipeline has been successfully completed:
1.  **Quality Control & Trimming:** Performed using `fastp` (Outputs in `fastp_qc/` and `multiqc_report/`).
2.  **rRNA Filtration:** Filtered out ribosomal RNA.
3.  **Mapping:** Reads were mapped to the reference genome (Outputs in `alignments/` and `filtered_alignments/`).
4.  **Quantification:** Absolute gene/transcript counts were generated using featureCounts (Located in `counts/gene_counts.txt`).

## 3. Downstream Analysis & Deliverables
The downstream analyses are found in the `downstream_results/` and `visualizations/` directories. 

### Data Tables (`downstream_results/`)
1.  **FPKM Normalized Counts** (`6_FPKM_normalized_counts_individual.csv`): Generated using DESeq2 without experimental design (`design = ~1`).
2.  **Annotated Expression Table** (`9_Annotated_Expression_Table_Individual.xlsx / .csv`): 
    *   The primary expression profiling table.
    *   Contains: Gene ID, Chromosome, Gene Name, Description, individual FPKM sample counts, Average Expression, Standard Deviation, and Z-scores.
    *   *Note: Gene Name and Description were extracted from `reference_genome/annotation.gtf`.*
3.  **Pairwise Comparisons** (`10_Pairwise_Comparisons.xlsx / .csv`):
    *   Contains Log2 Fold-Changes computed across all combinations of the 4 individual samples (e.g., NPFM1 vs NPFM2) on FPKM-normalized values. No p-values or FDR are included as there are no biological replicates.

### Visualizations (`visualizations/`)
All visualizations are generated natively in R/Python and exported as both `.pdf` and high-res `.png`.
1.  **PCA Plot** (`7_PCA_plot_individual`): Variance Stabilized Transformation (VST) based PCA. The background grid lines were explicitly removed via ggplot2 theme modifications.
2.  **Distance Matrix** (`8_Sample_Distance_Matrix_individual`): Euclidean distance clustering based on rlog counts.
3.  **Top 50 Variable Genes Heatmap** (`11_Heatmap_Top50_Variable_individual`): Expression signature of the 50 most variable genes across the 4 independent samples.
4.  **GO Classification Plot** (`GO_Classification_with_ontology_labels` & `without_ontology_labels`):
    *   Generated via `test_plotly_go.py`. 
    *   Y-axis is a log scale with the baseline lowered (`range=[-1.5, 2.05]`) so zero-value bars (mapped to `0.05`) are visually preserved.
    *   Extensive bottom margin (`b=500`) was added to give the ontology brackets (biological process, cellular component, molecular function) ample clearance (`y_offset=-1.0`).

## 4. GitHub & Web Deployment
*   **Repository:** The codebase, results, and scripts are synced to `https://github.com/Angsumi/MD-chakra`.
*   **Gitignore:** Massive raw `.fastq`, `.bam`, `.ht2`, and index files are explicitly ignored in `.gitignore` to adhere to GitHub file size limits.
*   **GitHub Pages (Web App):** 
    *   A premium frontend web page was built inside the `docs/` directory.
    *   **Live URL:** `https://angsumi.github.io/MD-chakra/`
    *   The page features glassmorphism design, serves embedded visualization PNGs, and provides functional download links (using absolute paths `/MD-chakra/assets/...`) for the core deliverables.
    *   Assets for the webpage are stored in `docs/assets/` and synced with the main results folders.

---
*Future Instruction for AGY: When resuming work on this repository, always read this file first to understand that no differential expression group-testing should be performed.*

## 5. Raw Data Access
*   **Google Drive Archive:** [Raw FASTQ, BAM, and Large Files](https://drive.google.com/drive/folders/1s6HS42fGcSoCgNkTawVGfZ5nhJUsGWVX?usp=drive_link)
    *   *Note: These files exceed GitHub's size limits and are hosted externally.*
