from rdkit import Chem
from rdkit.Chem import AllChem

class PeptideToSmilesConverter:
    def __init__(self, peptide_type="linear"):
        if peptide_type not in ["linear", "disulfide", "MBX", "KYL"]:
            raise ValueError("Invalid peptide_type. Choose 'linear', 'disulfide', 'MBX', or 'KYL'.")
        self.peptide_type = peptide_type

    def peptide_to_smiles(self, peptide_sequence):
        "Convert a peptide sequence to its corresponding SMILES representation."
        if self.peptide_type == "linear":
            return self._linear_peptide_to_smiles(peptide_sequence)
        elif self.peptide_type == "disulfide":
            return self._disulfide_peptide_to_smiles(peptide_sequence)
        elif self.peptide_type == "MBX":
            return self._mbx_peptide_to_smiles(peptide_sequence)

    def _linear_peptide_to_smiles(self, peptide_sequence):
        "Convert a linear peptide sequence to SMILES."
        mol = Chem.MolFromFASTA(peptide_sequence)
        if mol is None:
            raise ValueError(f"Unable to parse the peptide sequence: {peptide_sequence}")
        return Chem.MolToSmiles(mol)

    def _disulfide_peptide_to_smiles(self, peptide_sequence):
        mol = Chem.MolFromFASTA(peptide_sequence)
        if mol is None:
            raise ValueError(f"Unable to parse the peptide sequence: {peptide_sequence}")

        cysteine_indices = [i for i, aa in enumerate(peptide_sequence) if aa == 'C']
        if len(cysteine_indices) != 2:
            print(f"Skipping: {peptide_sequence} has {len(cysteine_indices)} cysteines (Expected exactly 2).")
            return None
        
        editable_mol = Chem.EditableMol(mol)
        sulfur_atoms = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'S':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetPDBResidueInfo() and neighbor.GetPDBResidueInfo().GetResidueName() == 'CYS':
                        sulfur_atoms.append(atom)
                        break
        if len(sulfur_atoms) != 2:
            print(f"Skipping: {peptide_sequence} has unexpected sulfur atom count: {len(sulfur_atoms)}")
            return None  # Return None to remove this row
            
        editable_mol.AddBond(sulfur_atoms[0].GetIdx(), sulfur_atoms[1].GetIdx(), Chem.BondType.SINGLE)
        mol_with_disulfide = editable_mol.GetMol()
        return Chem.MolToSmiles(mol_with_disulfide)

    def _mbx_peptide_to_smiles(self, peptide_sequence):
        "Convert linear peptide sequence to MBX-modified SMILES."
        linear_smiles = self._linear_peptide_to_smiles(peptide_sequence)
        mol = Chem.MolFromSmiles(linear_smiles)
        if mol is None:
            raise ValueError("Error converting sequence to SMILES.")

        modified_smiles = self._apply_mbx_modifications(mol, peptide_sequence)
        return modified_smiles

    def _apply_mbx_modifications(self, mol, peptide_sequence):
        """Modify peptides with exactly 2 cysteines (-SH) by adding an MBX (-CH₂-Benzene-CH₂-) bridge."""

        # Step 1: Identify cysteine residues using the peptide sequence
        cysteine_indices = [i for i, aa in enumerate(peptide_sequence) if aa == 'C']
    
        if len(cysteine_indices) != 2:
            print(f"Skipping MBX modification: {peptide_sequence} has {len(cysteine_indices)} cysteines.")
            return Chem.MolToSmiles(mol)

        # Step 2: Find sulfur atoms in the molecule
        sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']

        if len(sulfur_atoms) < 2:
            raise ValueError(f"Error: Expected at least 2 sulfur atoms, found {len(sulfur_atoms)} in {peptide_sequence}.")

        # Step 3: Find the hydrogen (H) atoms attached to sulfur (S) and remove them
        editable_mol = Chem.EditableMol(mol)
        sulfur_indices = [s.GetIdx() for s in sulfur_atoms]

        for sulfur in sulfur_atoms:
            for neighbor in sulfur.GetNeighbors():
                if neighbor.GetSymbol() == "H":  # Find the proton on sulfur
                    editable_mol.RemoveAtom(neighbor.GetIdx())

        # Step 4: Create the MBX linker (-CH₂-Benzene-CH₂-)
        mbx_fragment = Chem.MolFromSmiles("CC1=CC=CC=C1C")  # Example: Benzene (-CH₂-Ph-CH₂-)

    # Step 5: Attach the MBX linker to the two sulfur atoms
        combo = Chem.CombineMols(editable_mol.GetMol(), mbx_fragment)
        editable_combo = Chem.EditableMol(combo)

        s1_idx, s2_idx = sulfur_indices[0], sulfur_indices[1]
        benzene_c1 = mol.GetNumAtoms()  # First carbon in MBX
        benzene_c2 = mol.GetNumAtoms() + 7  # Second carbon in MBX

        editable_combo.AddBond(s1_idx, benzene_c1, Chem.BondType.SINGLE)
        editable_combo.AddBond(s2_idx, benzene_c2, Chem.BondType.SINGLE)

        new_mol = editable_combo.GetMol()
        return Chem.MolToSmiles(new_mol)
    
    def _apply_kyl_modifications(self, mol, peptide_sequence):
        """Modify peptides with exactly 2 cysteines (-SH) by adding the KYL (-S-CH2-Benzene-CH2-S-) bridge 
       and connecting an extra -CH2- group to the N-terminal amine of Serine."""
    
        # **Step 1: Identify cysteine residues using peptide sequence**
        cysteine_indices = [i for i, aa in enumerate(peptide_sequence) if aa == 'C']
    
        if len(cysteine_indices) != 2:
            print(f"Skipping KYL modification: {peptide_sequence} has {len(cysteine_indices)} cysteines.")
            return None  # ✅ Skip invalid sequences

        # **Step 2: Locate sulfur (S) atoms directly from cysteine positions**
        sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']
    
        if len(sulfur_atoms) < 2:
            print(f"Skipping KYL modification: {peptide_sequence} has fewer than 2 sulfur atoms.")
            return None  # ✅ Skip invalid sequences

        # **Step 3: Get sulfur atom indices for the two cysteines**
        s1_idx = sulfur_atoms[cysteine_indices[0]].GetIdx()
        s2_idx = sulfur_atoms[cysteine_indices[1]].GetIdx()

        print(f"Peptide: {peptide_sequence} - Connecting cysteine sulfurs at indices: {s1_idx}, {s2_idx}")

        # **Step 4: Find the N-terminal Serine Amine (`N`)**
        serine_index = peptide_sequence.find('S')  # Find Serine position
        if serine_index == -1:
            print(f"Skipping KYL modification: No serine found in {peptide_sequence}.")
            return None  # ✅ Skip if no serine

        # Find the amine (`N`) of the first Serine
        nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
        serine_n_idx = None
        for atom in nitrogen_atoms:
            if any(n.GetIdx() == serine_index for n in atom.GetNeighbors()):
                serine_n_idx = atom.GetIdx()
                break

        if serine_n_idx is None:
            print(f"Skipping KYL modification: Could not find Serine N-terminal amine in {peptide_sequence}.")
            return None  # ✅ Skip if no valid nitrogen

        print(f"Peptide: {peptide_sequence} - Serine N-terminal amine at index: {serine_n_idx}")

        # **Step 5: Create KYL bridge (-S-CH2-Benzene-CH2-S-) with extra -CH2- for Serine N-terminal amine**
        kyl_smiles = "CC1=CC(CN)=CC=C1CC"  # Extra -CH2- attached at ortho position
        kyl_fragment = Chem.MolFromSmiles(kyl_smiles)

        # **Step 6: Merge KYL bridge into original molecule**
        editable_mol = Chem.EditableMol(mol)
        combined_mol = Chem.CombineMols(editable_mol.GetMol(), kyl_fragment)
        combined_editable = Chem.EditableMol(combined_mol)

        # **Step 7: Identify KYL fragment attachment points**
        new_carbons = [atom for atom in combined_editable.GetMol().GetAtoms() if atom.GetSymbol() == "C"]

        if len(new_carbons) < 3:
            print(f"ERROR: Could not find expected KYL carbon atoms in {peptide_sequence}!")
            return None  # ✅ Skip if KYL not found

        kyl_c1_idx = new_carbons[-3].GetIdx()  # First benzyl carbon
        kyl_c2_idx = new_carbons[-2].GetIdx()  # Second benzyl carbon
        kyl_c3_idx = new_carbons[-1].GetIdx()  # Extra methylene (-CH2-)

        print(f"Peptide: {peptide_sequence} - KYL carbon atoms at indices: {kyl_c1_idx}, {kyl_c2_idx}, {kyl_c3_idx}")

        # **Step 8: Attach KYL linker to cysteine sulfur atoms**
        combined_editable.AddBond(s1_idx, kyl_c1_idx, Chem.BondType.SINGLE)
        combined_editable.AddBond(s2_idx, kyl_c2_idx, Chem.BondType.SINGLE)

    # **Step 9: Attach the extra methylene (-CH2-) to Serine N-terminal amine**
        combined_editable.AddBond(kyl_c3_idx, serine_n_idx, Chem.BondType.SINGLE)

        # **Step 10: Convert back to molecule**
        modified_mol = combined_editable.GetMol()
        modified_smiles = Chem.MolToSmiles(modified_mol)

        print(f"Peptide: {peptide_sequence} - Modified KYL SMILES: {modified_smiles}")

        return modified_smiles