import argparse
import os
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBList, Polypeptide
from Bio.SeqUtils.ProtParam import ProteinAnalysis

class Protein:
    """
    Encapsulates protein data and analysis methods.
    """
    def __init__(self, pdb_id: str, chain_id: str = None, pdb_dir: str = "pdb_files"):
        self.pdb_id = pdb_id
        self.chain_id = chain_id
        self.sequence = None
        self.structure = None
        self.pdb_dir = pdb_dir
        self._fetch_and_parse()

    def _fetch_and_parse(self):
        """Download PDB file and extract the amino acid sequence for the specified chain."""
        try:
            os.makedirs(self.pdb_dir, exist_ok=True)
            pdb_list = PDBList()
            pdb_file = pdb_list.retrieve_pdb_file(
                self.pdb_id, pdir=self.pdb_dir, file_format='pdb'
            )
            parser = PDBParser(QUIET=True)
            self.structure = parser.get_structure(self.pdb_id, pdb_file)
            chain_found = False
            for chain in self.structure[0]:
                if self.chain_id is None or chain.id == self.chain_id:
                    seq = Polypeptide.Polypeptide(chain).get_sequence()
                    self.sequence = str(seq)
                    chain_found = True
                    break
            if not chain_found or not self.sequence:
                raise ValueError("No valid chain/sequence extracted from PDB file.")
        except Exception as e:
            print(f"Error fetching/parsing PDB ID {self.pdb_id}: {e}")
            self.sequence = None

    def analyze_composition(self):
        """Print amino acid composition and plot as a bar chart."""
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
        plt.title(f"Amino Acid Composition ({self.pdb_id})")
        plt.xlabel("Amino Acid")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(f"{self.pdb_id}_composition.png")
        print(f"Saved amino acid composition bar chart to {self.pdb_id}_composition.png")

    def plot_hydrophobicity(self, window=9, show=False):
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
        plt.title(f"Hydrophobicity Plot for PDB ID: {self.pdb_id}")
        plt.xlabel("Residue Number")
        plt.ylabel("Kyte-Doolittle Hydrophobicity")
        plt.grid(True)
        plt.tight_layout()
        if show:
            plt.show()
        else:
            plt.savefig(f"{self.pdb_id}_hydrophobicity.png")
            print(f"Saved hydrophobicity plot to {self.pdb_id}_hydrophobicity.png")

    def get_info(self):
        """Return dictionary of basic protein properties."""
        if not self.sequence:
            return {"status": "Error", "message": "No sequence found."}
        analysis = ProteinAnalysis(self.sequence)
        info = {
            "PDB ID": self.pdb_id,
            "Chain": self.chain_id if self.chain_id else "First",
            "Sequence Length": len(self.sequence),
            "Molecular Weight": f"{analysis.molecular_weight():.2f} Da",
            "Instability Index": f"{analysis.instability_index():.2f}",
            "Aromaticity": f"{analysis.aromaticity():.2f}",
            "Isoelectric Point": f"{analysis.isoelectric_point():.2f}",
            "GRAVY": f"{analysis.gravy():.2f}"
        }
        return info

def main():
    parser = argparse.ArgumentParser(description="Protein Analysis Tool")
    parser.add_argument("pdb_id", help="Protein Data Bank ID (e.g., 1FAT)")
    parser.add_argument("--chain", type=str, default=None, help="Chain ID (e.g., A)")
    parser.add_argument("--window", type=int, default=9, help="Hydrophobicity window size")
    parser.add_argument("--show", action="store_true", help="Show plots interactively instead of saving")
    args = parser.parse_args()

    print(f"Analyzing PDB ID: {args.pdb_id}, Chain: {args.chain or 'First'}")
    protein = Protein(args.pdb_id, chain_id=args.chain)
    if not protein.sequence:
        print("Protein sequence extraction failed. Exiting.")
        return

    print("\n--- Protein Information ---")
    info = protein.get_info()
    for key, value in info.items():
        print(f"{key}: {value}")

    protein.analyze_composition()
    protein.plot_hydrophobicity(window=args.window, show=args.show)

if __name__ == "__main__":
    main()
