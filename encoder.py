from rdkit import Chem

class smilesEncoder:
    def __init__(self):
        self.visited = []

    def encode_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        self.visited = [False for _ in range(mol.GetNumAtoms())]
        newOrder =self. dfs(0, mol, [])
        nm = Chem.RenumberAtoms(mol, newOrder)
        return Chem.MolToSmiles(nm, canonical=False)

    def dfs(self, idx, mol, order):
        if self.visited[idx]:
            return order
        order.append(idx)
        self.visited[idx] = True
        atom = mol.GetAtomWithIdx(idx)
        nbr = [x.GetIdx() for x in atom.GetNeighbors()]
        #print(atom.GetSymbol(), atom.GetChiralTag(), nbr)
        hasFrom = False
        if nbr[0] < atom.GetIdx():
            hasFrom = True
        #print(hasFrom)
        if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
            if hasFrom:  # e.g. N[C@@](Br)(C)(O)
                nbr[1:] = reversed(nbr[1:])
            elif not hasFrom and len(nbr) == 4: # e.g. [C@H](N)(O)(C)
                nbr[1:] = reversed(nbr[1:])
        elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
            if not hasFrom and len(nbr) == 3:  # e.g. [C@H](Br)(C)O
                nbr = reversed(nbr)
        # print(nbr)
        for n in nbr:
            order = self.dfs(n, mol, order)
        return order

a = smilesEncoder()
print("Supported case:")
smiles = "N[C@](O)(Br)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "N[C@@](O)(Br)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "N[C@H](O)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "N[C@@H](O)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "[C@H](N)(O)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "[C@@H](N)(O)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "[C@@](Br)(N)(O)C"
print(smiles, '->', a.encode_smiles(smiles))
print("Unsupported case:")
smiles = "[C@]1(Br)(Cl)CCCC(F)C1"
print(smiles, '->', a.encode_smiles(smiles))
