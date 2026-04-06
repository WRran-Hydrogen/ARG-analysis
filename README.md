# SBL-virus

Code and analysis scripts for SBL virus research.

## Repository Structure

```
├── Caculate_distance.py          # Calculate pairwise distances
├── TranscrpticalAbundance_heatmap.R  # Transcriptional abundance heatmap
├── gggenome_domain.R             # GGGenome domain visualization
├── mg94_pairwise_plot_functions.r # mg94 pairwise plot functions
├── motif_enrichment.r            # Motif enrichment analysis
├── selection_pressure.r          # Selection pressure analysis
├── Molecular_docking/            # Molecular docking scripts
├── buble-world-map/              # Bubble world map visualization
└── sankey_graph/                 # Sankey graph visualization
```

## Dependencies

- R (>= 4.0)
- Python (>= 3.8)
- Required R packages: `ggplot2`, `reshape2`, `pheatmap`
- Required Python packages: `biopython`, `pandas`

## Usage

See individual script files for usage instructions.

## Data

Supporting data is organized into a **Figshare Collection** with 13 themed datasets across 4 parts:

| Dataset | Description |
|---------|-------------|
| Global Distribution Data | Geographic prevalence and location data |
| Motif Data | Motif node and edge relationships |
| Phylogeny Data | Phylogenetic trees (NCLDV, PF families) |
| Prophage Data | Prophage identification and InPhage data |
| DomainDistance Data | Pairwise domain distance calculations |
| OriginalFile Data | FASTA files, seed sequences, FIMO results |
| SelectionPressure Data | MG94 selection pressure analysis |
| MolecularDocking Data | PDB structures, PLIP analysis, docking results |
| Gephi Network Data | Network visualization files |
| Protein Structure Data (Part 1) | Predicted tertiary structures (PDB) |
| Protein Structure Data (Part 2) | Predicted tertiary structures (PDB) |
| Protein Structure Data (Part 3) | Predicted tertiary structures (PDB) |
| Protein Structure Data (Part 4) | Predicted tertiary structures (PDB) |

**Figshare Collection**: https://figshare.com/collections/_/8403348

Citation: xiong, yao (2026). SBL-virus Supporting Data. figshare. Collection. https://figshare.com/collections/_/8403348

## License

MIT License
