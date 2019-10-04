from rdkit import Chem

def encode_smiles(smiles, isDebug=False):
    mol = Chem.MolFromSmiles(smiles)
    visited = [False for _ in range(mol.GetNumAtoms())]
    def dfs(idx, order):
        if visited[idx]:
            return order
        order.append(idx)
        visited[idx] = True
        atom = mol.GetAtomWithIdx(idx)
        nbr = [x.GetIdx() for x in atom.GetNeighbors()]
        hasFrom = False
        if nbr[0] < atom.GetIdx():
            hasFrom = True

        if isDebug:
            print("Idx: {} ({}), isInRing: {}, hasFrom: {}".format(atom.GetIdx(), atom.GetSymbol(),  atom.IsInRing(), hasFrom),
                  atom.GetChiralTag(), nbr)

        if atom.IsInRing():
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                if hasFrom:  # e.g. 
                    if len(nbr) == 4:
                        nbr[1:] = reversed(nbr[1:])
                    else:
                        pass
                else:
                    if len(nbr) == 4:
                        pass
                    else:
                        pass
            elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                if hasFrom:  # e.g. 
                    if len(nbr) == 4:
                        nbr[1:] = reversed(nbr[1:])
                    else:
                        nbr[1:] = reversed(nbr[1:])
                else:
                    if len(nbr) == 4:  # e.g. [C@H](Br)(C)O
                        nbr[1:] = reversed(nbr[1:])
        else:
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                if hasFrom:  # e.g. N[C@@](Br)(C)(O)
                    nbr[1:] = reversed(nbr[1:])
                elif not hasFrom and len(nbr) == 4: # e.g. [C@H](N)(O)(C)
                    nbr[1:] = reversed(nbr[1:])
            elif atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    if not hasFrom and len(nbr) == 3:  # e.g. [C@H](Br)(C)O
                        nbr = reversed(nbr)
        if isDebug:
            print(nbr)
        for n in nbr:
            order = dfs(n, order)
        return order

    newOrder = dfs(0, [])
    if isDebug:
        print('new order:', newOrder)
    nm = Chem.RenumberAtoms(mol, newOrder)
    return Chem.MolToSmiles(nm, canonical=False)

def validate(smiles, isDebug=False):
    if isDebug:
        print('Validate:', smiles)
    stereo_smiles = encode_smiles(smiles, isDebug)
    assert(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == Chem.MolToSmiles(Chem.MolFromSmiles(stereo_smiles)))
    print('@@' not in stereo_smiles, smiles, '->', stereo_smiles)

validate('N[C@](O)(Br)C')
validate('N[C@@](O)(Br)C')
validate('N[C@H](O)C')
validate('N[C@@H](O)C')
validate('[C@H](N)(O)C')
validate('[C@@H](N)(O)C')
validate('[C@](Br)(N)(O)C')
validate('[C@@](Br)(N)(O)C')
validate('[C@]1(Br)(Cl)CCCC(F)C1')
validate('[C@@]1(Br)(Cl)CCCC(F)C1')
validate('N1C[C@H]2CC=CC[C@@H]2C1')
validate('C[C@@H]1CC(Nc2cncc(-c3nncn3C)c2)C[C@@H](C)C1')
validate('N#Cc1ccc(-c2ccc(O[C@@H](C(=O)N3CCCC3)c3ccccc3)cc2)cc1')
validate('CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C1')
validate('N#CC1=C(SCC(=O)Nc2cccc(Cl)c2)N=C([O-])[C@H](C#N)C12CCCCC2')
validate('CC[NH+](CC)[C@](C)(CC)[C@H](O)c1cscc1Br')
validate('COc1ccc(C(=O)N(C)[C@@H](C)C/C(N)=N/O)cc1O')
validate('C[C@@H]1CN(C(=O)c2cc(Br)cn2C)CC[C@H]1[NH3+]')
validate('CCOc1ccc(OCC)c([C@H]2C(C#N)=C(N)N(c3ccccc3C(F)(F)F)C3=C2C(=O)CCC3)c1')
validate('Cc1ccc2nc(S[C@H](C)C(=O)NC3CCC(C)CC3)n(C)c(=O)c2c1')
validate('Cc1ccccc1C(=O)N1CCC2(CC1)C[C@H](c1ccccc1)C(=O)N2C')
validate('CC(C)Cc1nc(SCC(=O)NC[C@@H]2CCCO2)c2c(=O)n(C)c(=O)n(C)c2n1')
validate('Cc1ccc(CNC(=O)c2ccccc2NC(=O)[C@@H]2CC(=O)N(c3ccc(C)cc3)C2)cc1')
validate('CC(C)[C@@H](Oc1cccc(Cl)c1)C(=O)N1CCC(n2cccn2)CC1')
validate('CCN(CC)C(=O)C[C@@H](C)[NH2+][C@H](C)c1cccc(F)c1F')
validate('O=C(NCCNC(=O)N1C[C@H]2CC=CC[C@@H]2C1)c1cccnc1')
validate('Cc1ccc(N2CC[C@@H](NS(=O)(=O)c3ccccc3C)C2=O)cc1C')
validate('CC[C@H](C)C[C@@H](C)NC(=O)N1CCN(CC(=O)NC2CC2)CC1')
validate('CC(=O)Nc1c2n(c3ccccc13)C[C@](C)(C(=O)NC1CCCCC1)N(C1CCCCC1)C2=O')
validate('N#Cc1ccncc1NC[C@@H]1C[C@@]12CCc1ccccc12')
validate('Cc1cccn2c(=O)c(C(=O)NC[C@H]3CCO[C@@H]3C(C)C)cnc12')
validate('C[C@@H](NC(=O)COC(=O)/C=C/c1ccc(Cl)cc1)c1ccccc1')
validate('CCOC(=O)c1nnc2ccccc2c1N1CC[C@@H]([NH+](CC)CC)C1')
validate('O=C(O[C@H]1CCOC1)C1(c2ccc(Cl)c(Cl)c2)CCC1')
validate('CCC[NH2+][C@@H]1COC[C@H]1C(=O)NCc1cscc1C')
validate('CC(=O)c1ccc(S(=O)(=O)N2CCCC[C@H]2C)cc1')
