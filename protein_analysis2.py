import argparse
import os
import sys
from typing import Optional, Dict
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBList, Polypeptide
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json

class Protein:
    def __init__(self, pdb_id: str, chain_id: Optional[str] = None, pdb_dir: str = "pdb_files", verbose: bool = True):
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.sequence = None
        self.structure = None
        self.pdb_dir = pdb_dir
        self.verbose = verbose
        self.pdb_file = None
        self._fetch_and_parse()

    def _fetch_and_parse(self):
        """Download PDB file (if needed) and extract sequence for specified chain."""
        try:
            os.makedirs(self.pdb_dir, exist_ok=True)
            filename = os.path.join(self.pdb_dir, f"pdb{self.pdb_id.lower()}.ent")
            if not os.path.exists(filename):
                if self.verbose:
                    print(f"Downloading PDB file for {self.pdb_id} ...")
                pdb_list = PDBList()
                self.pdb_file = pdb_list.retrieve_pdb_file(
                    self.pdb_id, pdir=self.pdb_dir, file_format='pdb'
                )
            else:
                self.pdb_file = filename
                if self.verbose:
                    print(f"Using cached PDB file {filename}")
            parser = PDBParser(QUIET=True)
            self.structure = parser.get_structure(self.pdb_id, self.pdb_file)
            # List available chains if chain_id not specified
            if self.chain_id is None:
                chains = [chain.id for chain in self.structure[0]]
                print(f"Available chains: {', '.join(chains)}")
                self.chain_id = chains[0]
                print(f"Using first chain: {self.chain_id}")
            chain_found = False
            for chain in self.structure[0]:
                if chain.id == self.chain_id:
                    seq = Polypeptide.Polypeptide(chain).get_sequence()
                    self.sequence = str(seq)
                    chain_found = True
                    break
            if not chain_found or not self.sequence:
                raise ValueError(f"Chain '{self.chain_id}' not found or no sequence extracted.")
        except Exception as e:
            print(f"Error fetching/parsing PDB ID {self.pdb_id}: {e}")
            self.sequence = None

    def analyze_composition(self, output_dir=".", save_csv=False):
        """Print and plot amino acid composition. Export to CSV if requested."""
        if not self.sequence:
            print("No sequence to analyze.")
            return
        analysis = ProteinAnalysis(self.sequence)
        composition = analysis.count_amino_acids()
        print("\n--- Amino Acid Composition ---")
        for aa, count in composition.items():
            print(f"{aa}: {count}")
        # Bar chart
        plt.figure(figsize=(8, 4))
        plt.bar(composition.keys(), composition.values(), color='skyblue')
        plt.title(f"Amino Acid Composition ({self.pdb_id}, chain {self.chain_id})")
        plt.xlabel("Amino Acid")
        plt.ylabel("Count")
        plt.tight_layout()
        png_path = os.path.join(output_dir, f"{self.pdb_id}_{self.chain_id}_composition.png")
        plt.savefig(png_path, dpi=300)
        print(f"Saved amino acid composition bar chart to {png_path}")
        if save_csv:
            csv_path = os.path.join(output_dir, f"{self.pdb_id}_{self.chain_id}_composition.csv")
            with open(csv_path, "w") as f:
                f.write("AminoAcid,Count\n")
                for aa, count in composition.items():
                    f.write(f"{aa},{count}\n")
            print(f"Saved composition data to {csv_path}")

    def plot_hydrophobicity(self, window=9, output_dir=".", show=False):
        """Create and save (or show) hydrophobicity plot."""
        if not self.sequence:
            print("No sequence for hydrophobicity plot.")
            return
        analysis = ProteinAnalysis(self.sequence)
        kd = {
            'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8, 'G': -0.4,
            'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
            'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8, 'T': -0.7, 'V': 4.2,
            'W': -0.9, 'Y': -1.3
        }
        hydrophobicity_values = analysis.protein_scale(kd, window, edge=1)
        plt.figure(figsize=(10, 6))
        plt.plot(hydrophobicity_values, color='coral')
        plt.title(f"Hydrophobicity Plot for PDB ID: {self.pdb_id} (chain {self.chain_id})")
        plt.xlabel("Residue Number")
        plt.ylabel("Kyte-Doolittle Hydrophobicity")
        plt.grid(True)
        plt.tight_layout()
        if show:
            plt.show()
        else:
            png_path = os.path.join(output_dir, f"{self.pdb_id}_{self.chain_id}_hydrophobicity.png")
            plt.savefig(png_path, dpi=300)
            print(f"Saved hydrophobicity plot to {png_path}")

    def get_info(self) -> Dict[str, str]:
        """Return dictionary of basic protein properties."""
        if not self.sequence:
            return {"status": "Error", "message": "No sequence found."}
        analysis = ProteinAnalysis(self.sequence)
        info = {
            "PDB ID": self.pdb_id,
            "Chain": self.chain_id,
            "Sequence Length": len(self.sequence),
            "Molecular Weight": f"{analysis.molecular_weight():.2f} Da",
            "Instability Index": f"{analysis.instability_index():.2f}",
            "Aromaticity": f"{analysis.aromaticity():.2f}",
            "Isoelectric Point": f"{analysis.isoelectric_point():.2f}",
            "GRAVY": f"{analysis.gravy():.2f}"
        }
        return info

    def export_info(self, output_dir=".", as_json=True):
        """Export protein info as JSON or CSV."""
        info = self.get_info()
        basename = f"{self.pdb_id}_{self.chain_id}_info"
        if as_json:
            json_path = os.path.join(output_dir, f"{basename}.json")
            with open(json_path, "w") as f:
                json.dump(info, f, indent=2)
            print(f"Saved protein info to {json_path}")
        else:
            csv_path = os.path.join(output_dir, f"{basename}.csv")
            with open(csv_path, "w") as f:
                for key, value in info.items():
                    f.write(f"{key},{value}\n")
            print(f"Saved protein info to {csv_path}")

def main():
    parser = argparse.ArgumentParser(description="Protein Analysis Tool (Enhanced)")
    parser.add_argument("pdb_id", help="Protein Data Bank ID (e.g., 1FAT)")
    parser.add_argument("--chain", type=str, default=None, help="Chain ID (e.g., A)")
    parser.add_argument("--window", type=int, default=9, help="Hydrophobicity window size")
    parser.add_argument("--show", action="store_true", help="Show plots interactively instead of saving")
    parser.add_argument("--output", type=str, default="output", help="Output directory for plots/data")
    parser.add_argument("--info-csv", action="store_true", help="Export protein info as CSV")
    parser.add_argument("--composition-csv", action="store_true", help="Export composition as CSV")
    parser.add_argument("--quiet", action="store_true", help="Suppress verbose output")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    protein = Protein(args.pdb_id, chain_id=args.chain, pdb_dir="pdb_files", verbose=not args.quiet)
    if not protein.sequence:
        print("Protein sequence extraction failed. Exiting.")
        sys.exit(1)

    print("\n--- Protein Information ---")
    info = protein.get_info()
    for key, value in info.items():
        print(f"{key}: {value}")

    protein.export_info(output_dir=args.output, as_json=not args.info_csv)
    protein.analyze_composition(output_dir=args.output, save_csv=args.composition_csv)
    protein.plot_hydrophobicity(window=args.window, output_dir=args.output, show=args.show)

if __name__ == "__main__":
    main()
