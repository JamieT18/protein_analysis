# Protein Analysis Tool

A robust Python tool for fetching, analyzing, and visualizing protein data from the Protein Data Bank (PDB).  
Given a PDB ID, the script downloads the structure, extracts the amino acid sequence (by chain), computes many properties, and generates publication-ready plots for amino acid composition and hydrophobicity.  
Results can be exported as CSV or JSON.

---

## Features

- **Fetch PDB file** by ID (cached locally for speed)
- **Extract amino acid sequence** (any specified chain; lists all available chains if unspecified)
- **Analyze composition and properties:**
  - Sequence length
  - Molecular weight
  - Instability index
  - Aromaticity
  - Isoelectric point
  - GRAVY score
- **Visualizations:**
  - Amino acid composition bar chart (`.png`, high-res)
  - Hydrophobicity plot (Kyte-Doolittle scale, sliding window, `.png`)
- **Flexible output options:**
  - Save plots to a custom directory
  - Show plots interactively (`--show`)
  - Export protein info and composition as CSV or JSON
- **User-friendly command-line interface:**
  - Lists available chains if not specified
  - Descriptive messages and error handling
  - Quiet mode (`--quiet`) for scripting
  - Caches downloads in `pdb_files` directory

---

## Setup

### 1. Clone or download the repository

```sh
git clone <your_repo_url>
cd <your_repo_dir>
```

### 2. Create a Python virtual environment

```sh
python -m venv venv
```

#### Activate on macOS/Linux:
```sh
source venv/bin/activate
```

#### Activate on Windows:
```sh
venv\Scripts\activate
```

### 3. Install requirements

```sh
pip install biopython matplotlib
```

---

## Usage

### Run the analysis script

```sh
python protein_analysis.py <PDB_ID> [--chain CHAIN] [--window WINDOW] [--show] [--output DIR] [--info-csv] [--composition-csv] [--quiet]
```

#### Arguments

- `<PDB_ID>`: PDB code (e.g., `1FAT`)
- `--chain CHAIN`: Specify chain ID (e.g., `A`). If omitted, all available chains are listed, and the first is used.
- `--window WINDOW`: Sliding window size for hydrophobicity plot (default: 9)
- `--show`: Display plots interactively instead of saving `.png` files
- `--output DIR`: Output directory for plots and exported data (default: `output`)
- `--info-csv`: Export protein info as CSV (default: JSON)
- `--composition-csv`: Export amino acid composition as CSV (in addition to plot)
- `--quiet`: Suppress verbose output

#### Examples

Analyze PDB ID **1FAT**, first chain, default options:
```sh
python protein_analysis.py 1FAT
```

Analyze PDB ID **2HYY**, chain B, window size 7, show plots interactively, output to `results/`:
```sh
python protein_analysis.py 2HYY --chain B --window 7 --show --output results
```

Export info and composition as CSV:
```sh
python protein_analysis.py 1FAT --info-csv --composition-csv
```

Suppress progress messages (quiet mode):
```sh
python protein_analysis.py 1FAT --quiet
```

---

## Output

- **Amino acid composition bar chart:** `output/PDBID_CHAIN_composition.png`
- **Hydrophobicity plot:** `output/PDBID_CHAIN_hydrophobicity.png`
- **Protein info:** `output/PDBID_CHAIN_info.json` (or `.csv` with `--info-csv`)
- **Amino acid composition table:** `output/PDBID_CHAIN_composition.csv` (with `--composition-csv`)
- **Protein info and composition:** printed to console

---

## Notes

- PDB files are cached in the `pdb_files` directory for faster re-analysis.
- An internet connection is required for new PDB downloads.
- If chain is unspecified, available chains are listed.
- For large proteins, hydrophobicity calculation may take longer.
- For more advanced analysis (secondary structure, mmCIF support), see future roadmap.

---

## Author

JamieT18

---

## Example Workflow

```sh
python protein_analysis.py 1FAT --chain A --window 7 --output analysis_results --show
```

This will:
- Download/cached the PDB structure
- Extract sequence for chain A
- Print properties and composition
- Save plots in `analysis_results` and show them interactively
- Export info as JSON (default)

---

## Future Enhancements

- Support for mmCIF files
- Secondary structure extraction
- Sequence logo visualization
- Integration with UniProt or Pfam
