# MD Chakra Transcriptomics

A transcriptomic analysis pipeline and profiling project for 4 individual samples (`NPFM1`, `NPFM2`, `NPFM3`, `NPFM4`). This project treats each sample independently without biological replicates, mapping reads and delivering comprehensive annotated expression tables.

## 🌟 Features
- **Independent Sample Profiling:** Accurate calculation of FPKM normalized counts and pairwise Log2 Fold-Changes without false grouped differential expression.
- **Annotated Expression Tables:** Gene profiling combining FPKM, average expression, standard deviations, Z-scores, and metadata extracted from reference GTF.
- **High-Quality Visualizations:** Includes PCA, sample distance matrices, heatmaps of the top 50 variable genes, and tailored GO classification plots.
- **Web Portal:** Includes a responsive web interface (glassmorphism design) to explore the analyses and visualizations.

## 🛠️ Tech Stack
- **Bioinformatics Tools:** fastp (QC), Bowtie2/Hisat (Mapping), featureCounts (Quantification)
- **Data Analysis & Visualization:** Python (Plotly), R (DESeq2, ggplot2)
- **Frontend Presentation:** HTML/CSS (GitHub Pages)

## 📊 Data & Results
- **Downstream Results:** Located in `downstream_results/`, containing final expression tables and pairwise comparison metrics.
- **Visualizations:** Stored in `visualizations/`, featuring all plotted metrics in PDF and PNG.
- **Raw Data:** All massive raw `.fastq`, `.bam`, and index files are hosted externally via [Google Drive](https://drive.google.com/drive/folders/1s6HS42fGcSoCgNkTawVGfZ5nhJUsGWVX?usp=drive_link) and omitted from git to preserve repository performance.

## 🌐 Live Web App
Explore the results and download core deliverables at the live portal:
[**MD Chakra Transcriptomics Portal**](https://angsumi.github.io/MD-chakra/)
