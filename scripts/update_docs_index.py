import re

html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Transcriptomic Analysis Results & 20-Analysis Suite | MD Chakra</title>
    <meta name="description" content="Comprehensive Transcriptomic expression profiling and 20-Analysis Suite for Nephila pilipes (NPFM1-4). Includes assembly QC, structural annotations, secretome, venom catalogue, and interactive deliverables.">
    <link rel="stylesheet" href="style.css">
    <style>
        .tier-header {
            margin-top: 3.5rem;
            margin-bottom: 1.5rem;
            padding-bottom: 0.8rem;
            border-bottom: 1px solid var(--border-hairline);
            display: flex;
            align-items: center;
            gap: 1rem;
        }
        .tier-badge {
            font-size: 0.75rem;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            padding: 0.35rem 0.85rem;
            border-radius: 9999px;
            font-weight: 700;
        }
        .tier-1-badge { background: rgba(99, 102, 241, 0.15); color: #818cf8; border: 1px solid rgba(99, 102, 241, 0.35); }
        .tier-2-badge { background: rgba(6, 182, 212, 0.15); color: #22d3ee; border: 1px solid rgba(6, 182, 212, 0.35); }
        .tier-3-badge { background: rgba(236, 72, 153, 0.15); color: #f472b6; border: 1px solid rgba(236, 72, 153, 0.35); }
        .quick-nav {
            display: flex;
            flex-wrap: wrap;
            gap: 0.75rem;
            margin: 2rem 0;
            padding: 1.25rem;
            background: var(--bg-card);
            border: 1px solid var(--border-hairline);
            border-radius: 12px;
        }
        .quick-nav a {
            color: var(--text-secondary);
            text-decoration: none;
            font-size: 0.9rem;
            font-weight: 600;
            padding: 0.4rem 0.85rem;
            border-radius: 8px;
            background: rgba(255, 255, 255, 0.04);
            transition: all 0.2s ease;
        }
        .quick-nav a:hover {
            color: #fff;
            background: var(--accent-indigo);
        }
    </style>
</head>
<body>

    <header>
        <h1>Transcriptomic Analysis Portal</h1>
        <p class="subtitle">Complete expression profiling and 20-Analysis Multi-Tier Transcriptome Suite for <i>Nephila pilipes</i> (<strong>NPFM1</strong>, <strong>NPFM2</strong>, <strong>NPFM3</strong>, <strong>NPFM4</strong>).</p>
        <div style="margin-top: 2rem; display: flex; justify-content: center; gap: 1rem; flex-wrap: wrap;">
            <a href="https://drive.google.com/drive/folders/1s6HS42fGcSoCgNkTawVGfZ5nhJUsGWVX?usp=drive_link" class="btn" target="_blank">Access Raw FASTQ & BAM on Google Drive</a>
            <a href="#suite-20" class="btn btn-secondary">Explore 20-Analysis Suite</a>
        </div>
    </header>

    <main>

        <!-- Quick Navigation -->
        <div class="quick-nav">
            <span style="color: var(--text-muted); font-weight: 700; align-self: center; margin-right: 0.5rem;">Quick Jump:</span>
            <a href="#matrix-table">Master Feasibility Table</a>
            <a href="#core-highlights">Core Deliverables</a>
            <a href="#tier-1">Tier 1: Assembly Quality (1-7)</a>
            <a href="#tier-2">Tier 2: Functional Annotation (8-14)</a>
            <a href="#tier-3">Tier 3: Spider Biology (15-20)</a>
        </div>

        <!-- Section 1: Deliverables Feasibility & Download Table -->
        <h2 id="matrix-table" class="section-title">Deliverables & Feasibility Matrix</h2>
        <p style="color: var(--text-muted); margin-bottom: 1.5rem;">
            Comprehensive status breakdown and direct download directory for all transcriptomic deliverables and functional profiling modules.
        </p>

        <div class="table-container">
            <table class="data-table">
                <thead>
                    <tr>
                        <th style="width: 50px;">#</th>
                        <th style="width: 260px;">Deliverable / Requirement</th>
                        <th style="width: 170px;">Feasibility Status</th>
                        <th>Methodology & Analytical Implementation</th>
                        <th style="width: 260px;">Download / Action</th>
                    </tr>
                </thead>
                <tbody>
                    <!-- Item 1 -->
                    <tr>
                        <td><strong>1</strong></td>
                        <td>
                            <strong>Raw & Processed FASTQ Files</strong>
                            <div class="item-desc">Paired-end reads before & after adapter trimming</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Trimmed with <code>fastp</code> (poly-G tail removal, adapter clipping, Q20/Q30 quality filtering) per sample.</td>
                        <td>
                            <a href="https://drive.google.com/drive/folders/1s6HS42fGcSoCgNkTawVGfZ5nhJUsGWVX?usp=drive_link" class="btn-sm btn-secondary" target="_blank">Google Drive (Raw)</a>
                        </td>
                    </tr>

                    <!-- Item 2 -->
                    <tr>
                        <td><strong>2</strong></td>
                        <td>
                            <strong>QC Reports & rRNA Filtration</strong>
                            <div class="item-desc">fastp summaries & ribosomal RNA depletion</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>SortMeRNA database indexing & Bowtie2 depletion (<code>--un-conc-gz</code>) with MultiQC report aggregation.</td>
                        <td>
                            <a href="/MD-chakra/assets/raw_qc_multiqc_report.html" class="btn-sm" download>MultiQC QC Report</a>
                        </td>
                    </tr>

                    <!-- Item 3 -->
                    <tr>
                        <td><strong>3</strong></td>
                        <td>
                            <strong>Filtered Mapping Statistics</strong>
                            <div class="item-desc">Total Reads, Mapped %, Unique %, Multi %, Unmapped %</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>HISAT2 genomic alignment metrics for rRNA-depleted mRNA reads across NPFM1, NPFM2, NPFM3, and NPFM4.</td>
                        <td>
                            <a href="/MD-chakra/assets/filtered_mapping_multiqc_report.html" class="btn-sm" download>MultiQC Alignment</a>
                            <a href="/MD-chakra/assets/Table_1_Summary_Stats.csv" class="btn-sm btn-secondary" download>Table 1 CSV</a>
                        </td>
                    </tr>

                    <!-- Item 4 -->
                    <tr>
                        <td><strong>4</strong></td>
                        <td>
                            <strong>Gene / Transcript Absolute Counts</strong>
                            <div class="item-desc">Expression profiling, variability, PCA, and normalization</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Quantified via <code>featureCounts</code> on genomic features; mapped to reference GTF locus IDs and gene symbols.</td>
                        <td>
                            <a href="/MD-chakra/assets/gene_counts_with_names.xlsx" class="btn-sm" download>Counts (.xlsx)</a>
                            <a href="/MD-chakra/assets/gene_counts.txt" class="btn-sm btn-secondary" download>Raw .txt</a>
                        </td>
                    </tr>

                    <!-- Item 5 -->
                    <tr>
                        <td><strong>5</strong></td>
                        <td>
                            <strong>FPKM Normalized Counts</strong>
                            <div class="item-desc">Depth and gene-length normalized expression</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Calculated using <code>DESeq2::fpkm()</code> with intercept design (<code>design = ~1</code>) and GTF feature lengths.</td>
                        <td>
                            <a href="/MD-chakra/assets/6_FPKM_normalized_counts_individual.csv" class="btn-sm" download>Download FPKM (.csv)</a>
                        </td>
                    </tr>

                    <!-- Item 6 -->
                    <tr>
                        <td><strong>6</strong></td>
                        <td>
                            <strong>Principal Component Analysis (PCA)</strong>
                            <div class="item-desc">Multivariate clustering across individual samples</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Computed using variance-stabilized <code>rlog</code> transformed counts showing sample separation and variance axes.</td>
                        <td>
                            <a href="/MD-chakra/assets/7_PCA_plot_individual.png" class="btn-sm" target="_blank">View PNG</a>
                            <a href="/MD-chakra/assets/7_PCA_plot_individual.pdf" class="btn-sm btn-secondary" download>PDF</a>
                        </td>
                    </tr>

                    <!-- Item 7 -->
                    <tr>
                        <td><strong>7</strong></td>
                        <td>
                            <strong>Sample Distance Matrix Heatmap</strong>
                            <div class="item-desc">Euclidean sample-to-sample distance clustering</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Hierarchical clustering on Euclidean distances of rlog counts across the 4 independent samples.</td>
                        <td>
                            <a href="/MD-chakra/assets/8_Sample_Distance_Matrix_individual.png" class="btn-sm" target="_blank">View PNG</a>
                            <a href="/MD-chakra/assets/8_Sample_Distance_Matrix_individual.pdf" class="btn-sm btn-secondary" download>PDF</a>
                        </td>
                    </tr>

                    <!-- Item 8 -->
                    <tr>
                        <td><strong>8</strong></td>
                        <td>
                            <strong>Differential Expression Analysis (DESeq2)</strong>
                            <div class="item-desc">Two-group hypothesis testing (Control vs Treatment)</div>
                        </td>
                        <td><span class="badge-danger">✖ Not Statistically Valid</span></td>
                        <td>Requires biological replicates ($n \ge 2$ per condition). Replaced with pairwise fold changes without false p-values.</td>
                        <td><span class="btn-sm btn-na">N/A (Bypassed)</span></td>
                    </tr>

                    <!-- Item 9 -->
                    <tr>
                        <td><strong>9</strong></td>
                        <td>
                            <strong>Master Annotated Expression Table</strong>
                            <div class="item-desc">Gene ID, Chr, Name, Description, S1..S4, Avg, StdDev, Z-scores</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Primary profiling deliverable: incorporates FPKM expression, mean, variance, and row Z-scores across all 4 samples.</td>
                        <td>
                            <a href="/MD-chakra/assets/9_Annotated_Expression_Table_Individual.xlsx" class="btn-sm" download>Master (.xlsx)</a>
                            <a href="/MD-chakra/assets/9_Annotated_Expression_Table_Individual.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                        </td>
                    </tr>

                    <!-- Item 10 -->
                    <tr>
                        <td><strong>10</strong></td>
                        <td>
                            <strong>Detailed Coordinate & Fold-Change Table</strong>
                            <div class="item-desc">Coords, counts, FC, Log2FC (without p-value/FDR)</div>
                        </td>
                        <td><span class="badge-warning">⚠ Partially Possible</span></td>
                        <td>All coordinates, counts, pairwise FC, and Log2FC are fully computed across all 6 pairwise sample combinations.</td>
                        <td>
                            <a href="/MD-chakra/assets/10_Pairwise_Comparisons.xlsx" class="btn-sm" download>Download Table (.xlsx)</a>
                        </td>
                    </tr>

                    <!-- Item 11 -->
                    <tr>
                        <td><strong>11</strong></td>
                        <td>
                            <strong>Top 50 Most Variable Genes Heatmap & Table</strong>
                            <div class="item-desc">Heatmap & annotated data table of highest variance genes</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Top 50 genes selected by row variance across normalized counts; standardized by row Z-scores with dedicated tables.</td>
                        <td>
                            <a href="/MD-chakra/assets/Top_50_Variable_Genes_Annotated.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                            <a href="/MD-chakra/assets/Top_50_Variable_Genes_Annotated.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                            <a href="/MD-chakra/assets/11_Heatmap_Top50_Variable_individual.png" class="btn-sm" target="_blank">View PNG</a>
                            <a href="/MD-chakra/assets/11_Heatmap_Top50_Variable_individual.pdf" class="btn-sm btn-secondary" download>PDF</a>
                        </td>
                    </tr>

                    <!-- Top Expressed Genes -->
                    <tr>
                        <td><strong>⭐</strong></td>
                        <td>
                            <strong>Top Highly Expressed Genes Tables</strong>
                            <div class="item-desc">Top 50 & Top 100 highest expressing genes ranked by mean FPKM</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Identifies dominant basal and structural transcripts across NPFM1–NPFM4 with complete gene annotation and expression metrics.</td>
                        <td>
                            <a href="/MD-chakra/assets/Top_50_Highly_Expressed_Genes.xlsx" class="btn-sm" download>Top 50 (.xlsx)</a>
                            <a href="/MD-chakra/assets/Top_50_Highly_Expressed_Genes.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                            <a href="/MD-chakra/assets/Top_100_Highly_Expressed_Genes.xlsx" class="btn-sm" download>Top 100 (.xlsx)</a>
                            <a href="/MD-chakra/assets/Top_100_Highly_Expressed_Genes.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                        </td>
                    </tr>

                    <!-- GO Functional Classification -->
                    <tr>
                        <td><strong>🧬</strong></td>
                        <td>
                            <strong>GO Functional Classification Dataset</strong>
                            <div class="item-desc">Figure 2 Gene Ontology breakdown (BP, CC, MF) with unigene counts & percentages</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Unigene functional annotation across Biological Process, Cellular Component, and Molecular Function categories.</td>
                        <td>
                            <a href="/MD-chakra/assets/Figure2_GO_Classification_Data.xlsx" class="btn-sm" download>Data (.xlsx)</a>
                            <a href="/MD-chakra/assets/Figure2_GO_Classification_Data.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                            <a href="/MD-chakra/assets/GO_Classification_with_ontology_labels.png" class="btn-sm btn-secondary" target="_blank">View Plot</a>
                        </td>
                    </tr>

                    <!-- KEGG Pathway Functional Annotation -->
                    <tr>
                        <td><strong>🗺️</strong></td>
                        <td>
                            <strong>KEGG Pathway Functional Annotation</strong>
                            <div class="item-desc">31,649 mapped unigenes, 770 KEGG pathways, KO numbers, EC enzymes, and BRITE functional hierarchy</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Diamond search against eggNOG orthology database with complete metabolic, cellular, and signaling pathway mappings.</td>
                        <td>
                            <a href="/MD-chakra/assets/KEGG_Pathway_Annotation_Master.xlsx" class="btn-sm" download>Master (.xlsx)</a>
                            <a href="/MD-chakra/assets/KEGG_Pathway_Annotation_Master.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                            <a href="/MD-chakra/assets/KEGG_Pathway_Classification_Stats.xlsx" class="btn-sm btn-secondary" download>Stats (.xlsx)</a>
                        </td>
                    </tr>

                    <!-- Protein Domain Functional Analysis -->
                    <tr>
                        <td><strong>🧩</strong></td>
                        <td>
                            <strong>Pfam Protein Domain Functional Analysis</strong>
                            <div class="item-desc">28,106 annotated transcripts across Toxins, Proteases, Ion Channels, Protein Binding, Enzymes & Extracellular Matrix</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Comprehensive Pfam domain identification detailing toxin signatures (ShKT, Kunitz), proteases (Trypsin, A17), ion channels, and binding scaffolds.</td>
                        <td>
                            <a href="/MD-chakra/assets/Protein_Domain_Annotation_Master.xlsx" class="btn-sm" download>Master (.xlsx)</a>
                            <a href="/MD-chakra/assets/Protein_Domain_Annotation_Master.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                            <a href="/MD-chakra/assets/Protein_Domain_Classification_Stats.xlsx" class="btn-sm btn-secondary" download>Stats (.xlsx)</a>
                        </td>
                    </tr>

                    
                    <!-- Digestion Functional Enzyme Repertoire -->
                    <tr style="background: rgba(16, 185, 129, 0.05);">
                        <td><strong>🍽️</strong></td>
                        <td>
                            <strong>Digestion: Functional Enzyme Repertoire & Clustered Heatmap</strong>
                            <div class="item-desc">1,952 curated digestive enzymes classified across 4 digestive processes and 28 enzyme families with abundance matrices</div>
                        </td>
                        <td><span class="badge-success">✔ Complete</span></td>
                        <td>Clustered heatmap of Enzyme Families × Functional Processes with transcript counts and abundance, plus gene-level replicate profiles.</td>
                        <td>
                            <a href="/MD-chakra/assets/01_Digestion_Enzyme_Master_Annotation.xlsx" class="btn-sm" download>Master (.xlsx)</a>
                            <a href="/MD-chakra/assets/01_Digestion_Functional_Enzyme_Heatmap.pdf" class="btn-sm btn-secondary" download>Heatmap (.pdf)</a>
                            <a href="/MD-chakra/assets/01_Digestion_Family_by_Process_Matrix.csv" class="btn-sm btn-secondary" download>Matrix (.csv)</a>
                        </td>
                    </tr>

                    <!-- Targeted Biological Processes Candidate Analysis -->
                    <tr>
                        <td><strong>🎯</strong></td>
                        <td>
                            <strong>Targeted Physiological Processes Candidate Genes</strong>
                            <div class="item-desc">Candidate genes, GO annotations, and KEGG pathways for Digestion, Detoxification, Sex Determination, Immunity, Silk Production & Toxin Biogenesis</div>
                        </td>
                        <td><span class="badge-success">✔ Fully Possible</span></td>
                        <td>Multi-sheet workbook detailing 10,647 candidate transcript instances annotated with exact Gene Symbols, Pfam domains, GO IDs, KEGG pathways, and expression profiles.</td>
                        <td>
                            <a href="/MD-chakra/assets/Targeted_Biological_Processes_Candidate_Genes.xlsx" class="btn-sm" download>Candidates (.xlsx)</a>
                            <a href="/MD-chakra/assets/Targeted_Biological_Processes_Candidate_Genes.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                            <a href="/MD-chakra/assets/Targeted_Processes_Summary_Stats.xlsx" class="btn-sm btn-secondary" download>Summary (.xlsx)</a>
                        </td>
                    </tr>

                </tbody>
            </table>
        </div>

        <!-- Section 2: Core Data Deliverables Cards -->
        <h2 id="core-highlights" class="section-title">Core Deliverable Highlights</h2>
        <div class="grid">

            <!-- Card: Digestion Repertoire -->
            <div class="card" style="border-color: rgba(16, 185, 129, 0.4);">
                <div class="card-header">
                    <h3 style="color: #10b981;">Digestive Enzyme Repertoire</h3>
                    <span class="badge" style="background: rgba(16, 185, 129, 0.2); color: #10b981;">.xlsx / .csv / .pdf</span>
                </div>
                <p>Curated master repertoire of 1,952 digestive enzymes across Protein, Lipid, Carbohydrate/Chitin, and Nucleic Acid digestion with exact replicate abundance.</p>
                <div class="process">
                    <strong>Deliverables:</strong>
                    Clustered Heatmap (Families × Processes), Gene-level Replicate Heatmap, Multi-sheet Master Excel, and Family Summary Matrices.
                </div>
                <div style="display: flex; gap: 0.5rem; margin-top: auto; flex-wrap: wrap;">
                    <a href="/MD-chakra/assets/01_Digestion_Enzyme_Master_Annotation.xlsx" class="btn" style="flex: 1;" download>Master Excel</a>
                    <a href="/MD-chakra/assets/01_Digestion_Functional_Enzyme_Heatmap.pdf" class="btn btn-secondary" style="flex: 1;" download>Heatmap PDF</a>
                </div>
            </div>
            
            <!-- Card 1 -->
            <div class="card">
                <div class="card-header">
                    <h3>Annotated Expression Table</h3>
                    <span class="badge">.xlsx / .csv</span>
                </div>
                <p>The master expression profiling matrix for all 4 samples. Contains Gene IDs, Chromosomal coordinates, Gene Names, Functional Descriptions, Average Expression, Standard Deviations, and Z-scores.</p>
                <div class="process">
                    <strong>Generation Process:</strong>
                    Quantified with featureCounts, normalized for sequencing depth and gene length, and merged with GTF attributes via pandas.
                </div>
                <div style="display: flex; gap: 0.5rem; margin-top: auto;">
                    <a href="/MD-chakra/assets/9_Annotated_Expression_Table_Individual.xlsx" class="btn" style="flex: 1;" download>Master (.xlsx)</a>
                    <a href="/MD-chakra/assets/9_Annotated_Expression_Table_Individual.csv" class="btn btn-secondary" style="flex: 1;" download>Master (.csv)</a>
                </div>
            </div>

            <!-- Card 2 -->
            <div class="card">
                <div class="card-header">
                    <h3>Pairwise Fold-Changes</h3>
                    <span class="badge">.xlsx</span>
                </div>
                <p>Pairwise Log2 Fold-Changes between all combinations of the 4 independent samples (e.g., NPFM1 vs NPFM2) to inspect relative expression changes without false group assumptions.</p>
                <div class="process">
                    <strong>Generation Process:</strong>
                    Direct logarithmic ratios calculated on FPKM-normalized values with pseudocounts across all 6 pairwise combinations.
                </div>
                <a href="/MD-chakra/assets/10_Pairwise_Comparisons.xlsx" class="btn" download>Download Pairwise Comparisons</a>
            </div>

            <!-- Card 3 -->
            <div class="card">
                <div class="card-header">
                    <h3>FPKM Normalized Counts</h3>
                    <span class="badge">.csv</span>
                </div>
                <p>Gene counts normalized for sequencing depth and gene length using Fragments Per Kilobase of transcript per Million mapped reads (FPKM).</p>
                <div class="process">
                    <strong>Generation Process:</strong>
                    Estimated size factors with DESeq2 treating samples as independent (<code>design = ~1</code>) factoring in feature lengths.
                </div>
                <a href="/MD-chakra/assets/6_FPKM_normalized_counts_individual.csv" class="btn" download>Download FPKM CSV</a>
            </div>

        </div>

        <!-- ============================================================= -->
        <!-- Section 3: Comprehensive 20-Analysis Transcriptome Portfolio -->
        <!-- ============================================================= -->
        <h2 id="suite-20" class="section-title">Comprehensive 20-Analysis Transcriptome Suite</h2>
        <p style="color: var(--text-muted); margin-bottom: 2rem;">
            Full systematic profiling covering Quality Control & Assembly Metrics (Tier 1), Functional & Structural Annotations (Tier 2), and Spider Biological Repertoires (Tier 3) for <i>Nephila pilipes</i>.
        </p>

        <!-- TIER 1 -->
        <div id="tier-1" class="tier-header">
            <span class="tier-badge tier-1-badge">Tier 1</span>
            <h3 style="font-size: 1.5rem; font-weight: 700; color: #f8fafc;">Assembly & Sequence Quality Statistics (Items 1 – 7)</h3>
        </div>
        
        <div class="grid">
            
            <!-- Analysis 1: BUSCO -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/01_BUSCO_Completeness.png" alt="BUSCO Completeness Assessment">
                </div>
                <div class="viz-content">
                    <h3>1. BUSCO Completeness</h3>
                    <p>Benchmarking Universal Single-Copy Orthologs assessment against Arthropoda and Arachnida lineages (>96% complete single & duplicated orthologs).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/01_BUSCO_Completeness.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/01_BUSCO_Completeness_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 2: Assembly Statistics -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/02_Assembly_Statistics_Lollipop.png" alt="Assembly Statistics Lollipop Chart">
                </div>
                <div class="viz-content">
                    <h3>2. Assembly Overview Statistics</h3>
                    <p>Lollipop chart profiling total transcript count (72,441), total length (82.1 Mb), mean/median lengths, and average GC content.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/02_Assembly_Statistics_Lollipop.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/02_Transcriptome_Assembly_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 3: N50/N90/L50/L90 -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/03_N50_N90_L50_L90_Metrics.png" alt="N50 N90 L50 L90 Contiguity Metrics">
                </div>
                <div class="viz-content">
                    <h3>3. N50/N90/L50/L90 Contiguity</h3>
                    <p>Contiguity metrics comparing contig length milestones (N50, N75, N90) and cumulative transcript count thresholds (L50, L75, L90).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/03_N50_N90_L50_L90_Metrics.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/03_Nx_Lx_Contiguity_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 4: Transcript Length ECDF -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/04_Transcript_Length_ECDF.png" alt="Transcript Length ECDF Plot">
                </div>
                <div class="viz-content">
                    <h3>4. Transcript Length Distribution</h3>
                    <p>Empirical Cumulative Distribution Function (ECDF) modeling transcript length dynamic range with median and N50 markers.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/04_Transcript_Length_ECDF.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 5: GC-Content Density -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/05_GC_Content_Density.png" alt="GC Content Density Profile">
                </div>
                <div class="viz-content">
                    <h3>5. GC-Content Distribution</h3>
                    <p>Kernel density estimation curve detailing GC percentage distribution across all 72,441 coding sequences (Mean: 41.2% GC).</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/05_GC_Content_Density.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 6: ORF / CDS Prediction -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/06_ORF_CDS_Prediction_Breakdown.png" alt="ORF CDS Prediction Breakdown">
                </div>
                <div class="viz-content">
                    <h3>6. ORF / CDS Completeness</h3>
                    <p>Classification of predicted open reading frames into Complete full-length CDS, 5'-partial, 3'-partial, and internal coding segments.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/06_ORF_CDS_Prediction_Breakdown.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/06_ORF_CDS_Prediction_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 7: Predicted Protein Stats -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/07_Predicted_Protein_Stats_ECDF.png" alt="Predicted Protein Stats ECDF">
                </div>
                <div class="viz-content">
                    <h3>7. Predicted Protein Statistics</h3>
                    <p>Cumulative length distribution (amino acids), molecular weights (kDa), and isoelectric point (pI) profiles across all translated proteins.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/07_Predicted_Protein_Stats_ECDF.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

        </div>

        <!-- TIER 2 -->
        <div id="tier-2" class="tier-header">
            <span class="tier-badge tier-2-badge">Tier 2</span>
            <h3 style="font-size: 1.5rem; font-weight: 700; color: #f8fafc;">Functional & Structural Annotation (Items 8 – 14)</h3>
        </div>

        <div class="grid">
            
            <!-- Analysis 8: InterPro Co-occurrence -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/08_InterPro_Domain_Cooccurrence.png" alt="InterPro Domain Co-occurrence Heatmap">
                </div>
                <div class="viz-content">
                    <h3>8. InterPro Domain Co-occurrence</h3>
                    <p>Heatmap demonstrating multi-domain architecture co-occurrence frequencies across top regulatory and structural domains.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/08_InterPro_Domain_Cooccurrence.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/08_InterPro_Domain_Cooccurrence_Matrix.xlsx" class="btn-sm" download>Matrix (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 9: eggNOG / COG Bubble Plot -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/09_eggNOG_COG_Classification_Bubble.png" alt="eggNOG COG Classification Bubble Plot">
                </div>
                <div class="viz-content">
                    <h3>9. eggNOG / COG Classification</h3>
                    <p>Bubble plot displaying unigene counts across all 25 standard COG categories (Information Storage, Cellular Processes, and Metabolism).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/09_eggNOG_COG_Classification_Bubble.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/09_eggNOG_COG_Functional_Classification.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 10: Transcription Factor Families -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/10_Transcription_Factor_Families.png" alt="Transcription Factor Families Lollipop Chart">
                </div>
                <div class="viz-content">
                    <h3>10. Transcription-Factor Families</h3>
                    <p>Ranked lollipop chart detailing all major TF families (C2H2-ZF, Homeobox, bHLH, bZIP, Forkhead, Nuclear Receptors, and DMRT).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/10_Transcription_Factor_Families.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/10_Transcription_Factor_Families_Stats.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 11: Signal Peptide Analysis -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/11_Signal_Peptide_Analysis_Donut.png" alt="Signal Peptide Analysis Donut Chart">
                </div>
                <div class="viz-content">
                    <h3>11. Signal Peptide Prediction</h3>
                    <p>SignalP donut chart showing proportion of proteins with cleavable N-terminal secretory signal peptides (SP+ vs SP-).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/11_Signal_Peptide_Analysis_Donut.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/11_Signal_Peptide_Prediction_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 12: Transmembrane Domains -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/12_Transmembrane_Domains_TMD.png" alt="Transmembrane Domains TMD Classes">
                </div>
                <div class="viz-content">
                    <h3>12. Transmembrane-Domain Analysis</h3>
                    <p>Distribution of TMHMM predicted transmembrane helices across classes (0 TMD soluble, 1 TMD single-pass, 7 TMD GPCRs, 8+ TMD channels).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/12_Transmembrane_Domains_TMD.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/12_Transmembrane_Domain_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 13: Secretome Prediction -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/13_Secretome_Prediction_UpSet.png" alt="Secretome Prediction Functional Classification">
                </div>
                <div class="viz-content">
                    <h3>13. Secretome Functional Repertoire</h3>
                    <p>Functional categorization of 8,920 high-confidence secreted proteins (SP+ / TM-) across toxins, proteases, structural silk, and effectors.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/13_Secretome_Prediction_UpSet.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/13_Secretome_Functional_Classification.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 14: Subcellular Localization -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/14_Subcellular_Localization.png" alt="Subcellular Localization Predictions">
                </div>
                <div class="viz-content">
                    <h3>14. Subcellular Localization</h3>
                    <p>DeepLoc compartmental localization predictions (Cytoplasm, Nucleus, Extracellular, Plasma Membrane, Mitochondrion, and ER).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/14_Subcellular_Localization.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/14_Subcellular_Localization_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

        </div>

        <!-- TIER 3 -->
        <div id="tier-3" class="tier-header">
            <span class="tier-badge tier-3-badge">Tier 3</span>
            <h3 style="font-size: 1.5rem; font-weight: 700; color: #f8fafc;">Spider Biology & Physiological Repertoires (Items 15 – 20)</h3>
        </div>

        <div class="grid">
            
            <!-- Analysis 15: Venom / Toxin Catalogue -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/15_Venom_Toxin_Catalogue.png" alt="Venom and Toxin Repertoire Catalogue">
                </div>
                <div class="viz-content">
                    <h3>15. Venom / Toxin Catalogue</h3>
                    <p>Repertoire of spider toxins and venom peptides: ShK potassium channel toxins, Kunitz inhibitors, Latrotoxin-like, and PLA2 enzymes.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/15_Venom_Toxin_Catalogue.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/15_Venom_Toxin_Repertoire_Catalogue.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 16: Protease Repertoire -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/16_Protease_Repertoire.png" alt="MEROPS Protease Repertoire">
                </div>
                <div class="viz-content">
                    <h3>16. Protease Repertoire</h3>
                    <p>Complete MEROPS protease family profiling: Aspartic (A17), Serine (S1 Trypsins), Cysteine (C1 Cathepsins, Caspases), and Astacin Metalloproteases.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/16_Protease_Repertoire.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/16_Protease_Repertoire_Catalogue.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 17: Ion-Channel Repertoire -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/17_Ion_Channel_Repertoire.png" alt="Ion Channel and Neuroreceptor Repertoire">
                </div>
                <div class="viz-content">
                    <h3>17. Ion-Channel Repertoire</h3>
                    <p>Voltage-gated (Kv, Cav, Nav) and ligand-gated (GluR, Nicotinic AChR, TRP, and ASIC/DEG acid-sensing) neuroreceptors and ion channels.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/17_Ion_Channel_Repertoire.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/17_Ion_Channel_Repertoire_Catalogue.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 18: Extracellular / Silk Repertoire -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/18_Extracellular_Structural_Silk.png" alt="Spidroin Silk Repertoire">
                </div>
                <div class="viz-content">
                    <h3>18. Spidroin Silk & Structural Matrix</h3>
                    <p>Full spider silk spidroin catalog: MaSp1/2 (dragline), MiSp (web frame), Flag (capture spiral), TuSp (egg sac), PySp, and cuticle chitin-binding proteins.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/18_Extracellular_Structural_Silk.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/18_Silk_and_Structural_Protein_Catalogue.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 19: Reference Genes Stability -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/19_Candidate_Reference_Genes_Stability.png" alt="Candidate Reference Genes Stability Heatmap">
                </div>
                <div class="viz-content">
                    <h3>19. Candidate Reference Genes</h3>
                    <p>geNorm stability score (M) and coefficient of variation (CV%) ranking for internal RT-qPCR reference gene selection (GAPDH, RPL32, EF1A, ACTB).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/19_Candidate_Reference_Genes_Stability.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/19_Candidate_Reference_Genes_Stability.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Analysis 20: Alternative Splicing -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/20_Alternative_Splicing_Isoforms.png" alt="Alternative Splicing Events Profile">
                </div>
                <div class="viz-content">
                    <h3>20. Alternative Splicing / Isoforms</h3>
                    <p>Genome-wide alternative splicing profiling: Alternative 3' splice sites (A3SS), Alternative 5' splice sites (A5SS), Skipped Exons (SE), and Retained Introns (RI).</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <a href="/MD-chakra/assets/20_Alternative_Splicing_Isoforms.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/20_Alternative_Splicing_Events_Statistics.xlsx" class="btn-sm" download>Table (.xlsx)</a>
                    </div>
                </div>
            </div>

        </div>

        <!-- ============================================================= -->
        <!-- Section 4: Multivariate Expression Visualizations -->
        <!-- ============================================================= -->
        <h2 class="section-title" style="margin-top: 4rem;">Multivariate Expression & Functional Visualizations</h2>
        <div class="grid">

            <!-- Viz: Digestion Clustered Heatmap -->
            <div class="viz-card" style="border-color: rgba(16, 185, 129, 0.35);">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/01_Digestion_Functional_Enzyme_Heatmap.png" alt="Digestion Functional Enzyme Clustered Heatmap">
                </div>
                <div class="viz-content">
                    <h3>Functional Digestive Enzyme Heatmap</h3>
                    <p>Hierarchically clustered heatmap displaying Enzyme Families × Functional Digestive Processes with transcript count overlays and log2-abundance shading.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem;">
                        <a href="/MD-chakra/assets/01_Digestion_Functional_Enzyme_Heatmap.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/01_Digestion_Functional_Enzyme_Heatmap.png" class="btn-sm btn-secondary" download>Download PNG</a>
                    </div>
                </div>
            </div>

            <!-- Viz: Digestion Gene Expression Heatmap -->
            <div class="viz-card" style="border-color: rgba(59, 130, 246, 0.35);">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/01b_Digestion_Gene_Expression_Heatmap.png" alt="Top Expressed Digestive Enzymes Across Replicates">
                </div>
                <div class="viz-content">
                    <h3>Top Expressed Digestive Enzymes</h3>
                    <p>Sample-level clustered expression profiles (NPFM1–NPFM4) of top candidate digestive enzymes with enzyme class sidebar color-coding.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem;">
                        <a href="/MD-chakra/assets/01b_Digestion_Gene_Expression_Heatmap.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/01b_Digestion_Gene_Expression_Heatmap.png" class="btn-sm btn-secondary" download>Download PNG</a>
                    </div>
                </div>
            </div>

            <!-- Viz: Digestion Repertoire Breakdown -->
            <div class="viz-card" style="border-color: rgba(245, 158, 11, 0.35);">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/01c_Digestion_Family_Repertoire_Breakdown.png" alt="Digestion Enzyme Family Diversity & Abundance Profile">
                </div>
                <div class="viz-content">
                    <h3>Digestive Enzyme Diversity & Abundance</h3>
                    <p>Comparative dual-panel profile of gene repertoire richness versus total transcriptional output across digestive functional classes.</p>
                    <div style="margin-top: 1rem; display: flex; gap: 0.5rem;">
                        <a href="/MD-chakra/assets/01c_Digestion_Family_Repertoire_Breakdown.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/01c_Digestion_Family_Repertoire_Breakdown.png" class="btn-sm btn-secondary" download>Download PNG</a>
                    </div>
                </div>
            </div>

            <!-- Viz 1: PCA -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/7_PCA_plot_individual.png" alt="Principal Component Analysis Plot">
                </div>
                <div class="viz-content">
                    <h3>PCA Plot</h3>
                    <p>Principal Component Analysis based on Variance Stabilized Transformation (rlog) of read counts across individual samples.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/7_PCA_plot_individual.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 2: Distance Matrix -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/8_Sample_Distance_Matrix_individual.png" alt="Sample Distance Matrix">
                </div>
                <div class="viz-content">
                    <h3>Sample Distance Matrix</h3>
                    <p>Euclidean distance clustering of the 4 samples based on rlog-transformed counts.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/8_Sample_Distance_Matrix_individual.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 3: Top 50 Variable -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/11_Heatmap_Top50_Variable_individual.png" alt="Top 50 Variable Genes Heatmap">
                </div>
                <div class="viz-content">
                    <h3>Top 50 Variable Genes</h3>
                    <p>Heatmap illustrating the expression patterns of the 50 genes with the highest variance across samples, scaled by row Z-scores.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/11_Heatmap_Top50_Variable_individual.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>
            
            <!-- Viz 4: GO Classification -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/GO_Classification_with_ontology_labels.png" alt="GO Classification Plot">
                </div>
                <div class="viz-content">
                    <h3>GO Classification</h3>
                    <p>Gene Ontology (GO) classification bar chart displaying unigene counts across Biological Process, Cellular Component, and Molecular Function.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/Figure2_GO_Classification_Data.xlsx" class="btn-sm btn-secondary" download>Dataset (.xlsx)</a>
                        <a href="/MD-chakra/assets/Figure2_GO_Classification_Data.csv" class="btn-sm btn-secondary" download>(.csv)</a>
                    </div>
                </div>
            </div>

            <!-- Viz 5: Classic Tukey Boxplot -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/12_Expression_Distribution_Boxplot.png" alt="Log2-Expression Boxplot">
                </div>
                <div class="viz-content">
                    <h3>Log2-Expression Boxplot</h3>
                    <p>Classic Tukey boxplot displaying median (solid line), mean (dashed line), 25th–75th IQR, and whiskers of log2(FPKM + 1) across NPFM1–NPFM4.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/12_Expression_Distribution_Boxplot.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 5b: Standalone Violin Plot -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/12b_Expression_Distribution_Violin.png" alt="Log2-Expression Violin Distributions">
                </div>
                <div class="viz-content">
                    <h3>Log2-Expression Violin Plot</h3>
                    <p>Continuous violin density profiles demonstrating matching multi-modal distribution dynamics across all 4 libraries.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/12b_Expression_Distribution_Violin.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 6: Density Plot (KDE) -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/13_Expression_Density_Plot.png" alt="Expression Density Distribution KDE">
                </div>
                <div class="viz-content">
                    <h3>Expression Density Distributions (KDE)</h3>
                    <p>Kernel Density Estimation curves displaying matching multi-modal expression profiles across all 4 independent libraries.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/13_Expression_Density_Plot.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 7: Cumulative CPM Distribution -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/14_Cumulative_CPM_Distribution.png" alt="Cumulative CPM Read Fraction">
                </div>
                <div class="viz-content">
                    <h3>Cumulative CPM Read Fraction</h3>
                    <p>Quantifies library complexity and transcript diversity, tracking the cumulative proportion of sequencing depth consumed across unigenes.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/14_Cumulative_CPM_Distribution.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 8: KEGG Pathway Classification Bar Chart -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/15_KEGG_Pathway_Classification.png" alt="KEGG Pathway Functional Classification">
                </div>
                <div class="viz-content">
                    <h3>KEGG Pathway Classification</h3>
                    <p>Top 20 represented biological and signaling pathways across the spider transcriptome, highlighting Endocytosis, Ribosome, Lysosome, and Spliceosome.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/15_KEGG_Pathway_Classification.pdf" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/KEGG_Pathway_Annotation_Master.xlsx" class="btn-sm" download>Master Table (.xlsx)</a>
                    </div>
                </div>
            </div>

            <!-- Viz 9: KEGG Pathway Representation Dotplot -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/16_KEGG_Pathway_Dotplot.png?v=20260911_02" alt="KEGG Pathway Representation Dotplot">
                </div>
                <div class="viz-content">
                    <h3>KEGG Pathway Representation Dotplot</h3>
                    <p>Publication ggplot2-style dotplot illustrating pathway representation, gene ratios, and significance values.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/16_KEGG_Pathway_Dotplot.pdf?v=20260911_02" class="btn-sm btn-secondary" download>Download PDF</a>
                    </div>
                </div>
            </div>

            <!-- Viz 10: Protein Domain Functional Analysis -->
            <div class="viz-card">
                <div class="viz-img-container">
                    <img src="/MD-chakra/assets/17_Protein_Domain_Classification.png?v=20260911_01" alt="Top Protein Domains in Nephila pilipes Control Transcriptome">
                </div>
                <div class="viz-content">
                    <h3>Top Protein Domains Analysis</h3>
                    <p>6-panel categorized Pfam protein domain architecture highlighting Toxins, Proteases, Ion Channels, Binding, Enzymes, and Extracellular Matrix.</p>
                    <div style="margin-top: 1rem;">
                        <a href="/MD-chakra/assets/17_Protein_Domain_Classification.pdf?v=20260911_01" class="btn-sm btn-secondary" download>Download PDF</a>
                        <a href="/MD-chakra/assets/Protein_Domain_Annotation_Master.xlsx" class="btn-sm" download>Master Table (.xlsx)</a>
                    </div>
                </div>
            </div>

        </div>
    </main>

    <footer>
        <p>MD Chakra Transcriptomics Pipeline &copy; 2026.</p>
    </footer>

</body>
</html>
"""

with open('docs/index.html', 'w') as f:
    f.write(html_content)

print("Updated docs/index.html with full 20-Analysis Suite and interactive visuals!")
