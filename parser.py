from rdkit import Chem
import copy

def smilesToMol(smiles):
    prev = -1
    vprev = []
    rclose = []
    mol = Molecule()
    ptr = 0
    order = 0
    updown = ''

    def parseSimple():
        nonlocal mol, ptr, prev, vprev, order, updown, rclose
        element = smiles[ptr]
        if ptr+1 < len(smiles):
            if element == 'B' and smiles[ptr+1] == 'r':
                element += smiles[ptr+1]
                ptr += 1
            if element == 'C' and smiles[ptr+1] == 'l':
                element += smiles[ptr+1]
                ptr += 1

        atom = mol.AddAtom()
        atom.element = element
        if prev != -1:
            if order == 0:
                mol.AddBond(prev, atom.idx, 1)
            else:
                mol.AddBond(prev, atom.idx, order)
        prev = atom.idx
        order = 0
        return True

    def parseComplex():
        nonlocal mol, ptr, prev, vprev, order, updown, rclose
        ptr += 1
        atom = mol.AddAtom()
        element = smiles[ptr]
        if element.isupper() and ptr+1 < len(smiles) and smiles[ptr+1].islower():
            element += smiles[ptr+1]
            ptr += 1
        atom.element = element
        ptr += 1
        charge = 0
        hcount = 0
        while smiles[ptr] != ']':
            token = smiles[ptr]
            if token == '@':
                if smiles[ptr+1] == '@':
                    atom.tetrahedralStereo = 'Clockwise'
                    ptr += 1
                else:
                    atom.tetrahedralStereo = 'CounterClockwise'
            elif token == '-':
                ptr += 1
                if not smiles[ptr].isdigit():
                    charge -= 1
                while smiles[ptr].isdigit():
                    charge = charge*10 - int(smiles[ptr])
                    ptr += 1
                ptr -= 1
            elif token == '+':
                ptr += 1
                if not smiles[ptr].isdigit():
                    charge += 1
                while smiles[ptr].isdigit():
                    charge = charge*10 + int(smiles[ptr])
                    ptr += 1
                ptr -= 1
            elif token == 'H':
                ptr += 1
                if smiles[ptr].isdigit():
                    hcount = int(smiles[ptr])
                else:
                    hcount = 1
                    ptr -= 1
            else:
                return False
            atom.charge = charge
            atom.hcount = hcount
            ptr += 1

        if prev != -1:
            if order == 0:
                mol.AddBond(prev, atom.idx, 1)
            else:
                mol.AddBond(prev, atom.idx, order)
        prev = atom.idx
        return True

    def parseRingBond():
        nonlocal ptr, prev, vprev, order, updown, rclose
        prevAtom = mol.GetAtom(prev)
        if prevAtom is None:
            return False
        digit = int(smiles[ptr])
        for rc in rclose:
            if rc.digit == digit:
                # Check for self-bonding, e.g. C11
                if rc.prev == prev:
                    return False
                bondOrder = order if order > rc.order else rc.order
                if bondOrder == 0:
                    bondOrder = 1
                mol.AddBond(rc.prev, prev, bondOrder)
                rclose.remove(rc)
                updown = ''
                order = 0
                return True
        # no closure
        rclose.append(ringClosure(digit, prev, order))
        order = 0
        updown = ''
        return True

    def parseSmiles(smiles):
        nonlocal ptr, prev, vprev, order, updown
        while ptr < len(smiles):
            token = smiles[ptr]
            # ring open/close
            if token in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                if prev == -1:
                    return False
                if not parseRingBond():
                    return False
            elif token == '(':
                vprev.append(prev)
            elif token == ')':
                if not vprev:  # vprev is empty
                    return False
                prev = vprev.pop()
            elif token == '[':
                if not parseComplex():
                    return False
            elif token in ['-', '=', '#', '/', '\\']:
                if prev == -1:
                    return False
                if token == '-':
                    order = 1
                elif token == '=':
                    order = 2
                elif token == '#':
                    order = 3
                elif token == '/':
                    order = 1
                    updown = 'Down'
                elif token == '\\':
                    order = 1
                    updown = 'Up'
            else:
                if not parseSimple():
                    return False
            ptr += 1
        return True

    parseSmiles(smiles)
    return mol


class Molecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []

    def AddAtom(self):
        atom = Atom()
        atom.idx = len(self.atoms)
        self.atoms.append(atom)
        return atom

    def GetAtom(self, idx):
        if idx >= 0 and idx < len(self.atoms):
            return self.atoms[idx]
        else:
            return None

    def NumAtoms(self):
        return len(self.atoms)

    def AddBond(self, begin, end, order):
        bond = Bond(begin, end, order)
        self.bonds.append(bond)
        self.atoms[begin].nbr.append(end)
        self.atoms[end].nbr.append(begin)
        return bond

    def GetBond(self, idx):
        if idx >= 0 and idx < len(self.bonds):
            return self.bonds[idx]
        else:
            return None

    def NumBonds(self):
        return len(self.bonds)


class Atom:
    def __init__(self):
        self.idx = -1
        self.element = ""
        self.charge = 0
        self.hcount = 0
        self.newidx = -1
        self.nbr = []
        self.tetrahedralStereo = None
        self.ringClosure = False

    def __repr__(self):
        return "({}, idx: {}, newidx: {}, charge: {}, h: {}, nbr: {}, rc: {}, {})".format(self.element, self.idx, self.newidx, self.charge, self.hcount, self.nbr, self.ringClosure, self.tetrahedralStereo)


class Bond:
    def __init__(self, begin=-1, end=-1, order=1):
        self.begin = begin
        self.end = end
        self.order = order

    def __repr__(self):
        return "(begin: {}, end: {}, order: {})".format(
                self.begin, self.end, self.order)


class ringClosure:
    def __init__(self, digit, prev, order):
        self.digit = digit
        self.prev = prev
        self.order = order

    def __repr__(self):
        return "(digit: {}, prev: {}, order: {})".format(
                self.digit, self.prev, self.order)

# find ring closure
def generateSmiles(mol):
    # decide atom order and find ring closure
    order = []
    visited = [False for _ in range(len(mol.atoms))]
    stack = [(-1, 0)]
    rcbond = []

    while stack:
        prev, idx = stack.pop()
        if not visited[idx]:
            order.append(idx)
            mol.atoms[idx].newidx = len(order) - 1
            visited[idx] = True
            atom = mol.atoms[idx]
            if atom.tetrahedralStereo == 'Clockwise':
                if atom.nbr[2] == prev:
                    atom.nbr[:2] = [atom.nbr[1], atom.nbr[0]]
                else:
                    atom.nbr[-2:] = [atom.nbr[-1], atom.nbr[-2]]
            elif atom.tetrahedralStereo == 'CounterClockwise':
                if atom.nbr[1] == prev:
                    atom.nbr[0], atom.nbr[2] = atom.nbr[2], atom.nbr[0]
                
            nbr = reversed(atom.nbr)
            for a in nbr:
                if not visited[a]:
                    stack.append((idx, a))
                elif a != prev:
                    mol.atoms[a].ringClosure = True
                    rcbond.append((a, idx))

    usedDigit = []
    rcstack = []
    pstack = []
    #for atom in mol.atoms:
    #    print(atom)

    def generateBranchSmiles(idx):
        nonlocal mol, order, usedDigit, rcstack, rcbond, pstack
        smiles = ""
        prev = pstack.pop()
        while True:
            atom = mol.atoms[idx]
            if atom.tetrahedralStereo is None:
                if atom.charge == 0:
                    smiles += atom.element
                else:
                    if atom.charge == -1:
                        smiles += '[{}-]'.format(atom.element)
                    elif atom.charge == +1:
                        if atom.hcount == +1:
                            smiles += '[{}H+]'.format(atom.element)
                        elif atom.hcount == +2:
                            smiles += '[{}H2+]'.format(atom.element)
                        elif atom.hcount == +3:
                            smiles += '[{}H3+]'.format(atom.element)
                        else:
                            print('Unsupported!')
                            exit()
                    else:
                        print('Unsupported!')
                        exit()
            else:
                smiles += '[' + atom.element + '@'
                if atom.hcount == 1:
                    smiles += 'H'
                smiles += ']'
            # if ring is starting
            if atom.ringClosure:
                if not usedDigit:
                    digit = 1
                else:
                    digit = max(usedDigit) + 1
                usedDigit.append(digit)
                rc = ringClosure(digit, idx, -1)
                rcstack.append(rc)
                smiles += str(digit)
            # if there is branch (buggy)
            branch = [n for n in atom.nbr if n != prev and (idx, n) not in rcbond]
            if not branch:
                return idx, smiles
            while len(branch) > 1:
                bidx = branch.pop(0)
                # check ring closure
                isRingClosure = False
                for rc in rcstack:
                    if rc.prev == bidx:
                        smiles += str(rc.digit)
                        rcstack.remove(rc)
                        isRingClosure = True
                        break
                if not isRingClosure:
                    smiles += '('
                    for b in mol.bonds:
                        if b.begin == idx and b.end == bidx:
                            if b.order == 2:
                                smiles += '='
                            elif b.order == 3:
                                smiles += '#'
                    pstack.append(idx)
                    _, branch_smiles = generateBranchSmiles(bidx)
                    smiles += branch_smiles
                    smiles += ')'
            nidx = branch.pop()
            # if end of ring
            if mol.atoms[nidx].newidx < atom.newidx:
                if not rcstack:
                    return idx, smiles
                for rc in rcstack:
                    if rc.prev == nidx:
                        smiles += str(rc.digit)
                        rcstack.remove(rc)
                        usedDigit.remove(rc.digit)
                        return idx, smiles
                print("something is wrong!")
                exit()
            # ordinary case
            else:
                for b in mol.bonds:
                    if b.begin == idx and b.end == nidx:
                        if b.order == 2:
                            smiles += '='
                        elif b.order == 3:
                            smiles += '#'

                prev = idx
                idx = nidx
    pstack = [-1]
    _, smiles = generateBranchSmiles(0)
    return smiles

def validate(smiles):
    mol = smilesToMol(smiles)
    new_smiles = generateSmiles(mol)
    isCorrect = (Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == Chem.MolToSmiles(Chem.MolFromSmiles(new_smiles)))
    print(isCorrect, smiles, '->', new_smiles)

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
validate('CCOC(=O)c1nnc2ccccc2c1N1CC[C@@H]([NH+](CC)CC)C1')
validate('CCC[NH2+][C@@H]1COC[C@H]1C(=O)NCc1cscc1C')
validate('CC(=O)c1ccc(S(=O)(=O)N2CCCC[C@H]2C)cc1')
