from rdkit import Chem
import copy

def smilesToMol(smiles):
    # normalize
    smiles = smiles.replace('[C@H]', '[C@]([H])')
    smiles = smiles.replace('[C@@H]', '[C@@]([H])')

    prev = -1
    vprev = []
    rclose = []
    mol = Molecule()
    ptr = 0
    order = 0
    updown = ''

    def parseSimple():
        nonlocal mol, ptr, prev, vprev, order, updown, rclose, updown
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
        if updown == '/' or updown == '\\':
            mol.upDownMap[len(mol.bonds)-1] = updown
        updown = ''
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

        if updown == '/' or updown == '\\':
            mol.upDownMap[len(mol.bonds)-1] = updown
        updown = ''

        prev = atom.idx
        order = 0
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
                i = mol.atoms[rc.prev].nbr2.index(-digit)
                mol.atoms[rc.prev].nbr2[i] = mol.atoms[rc.prev].nbr2.pop()
                if updown == '/' or updown == '\\':
                    mol.upDownMap[len(mol.bonds)-1] = updown
                if rc.updown == '/' or rc.updown == '\\':
                    if len(mol.bonds) - 1 in mol.upDownMap:
                        print("Illegal cis / trans specification")
                        exit()
                    mol.upDownMap[len(mol.bonds)-1] = updown
                updown = ''
                rclose.remove(rc)
                order = 0
                return True
        # no closure
        rclose.append(ringClosure(digit, prev, order, updown))
        mol.atoms[prev].nbr2.append(-digit)
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
                    updown = '/'
                elif token == '\\':
                    order = 1
                    updown = '\\'
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
        self.upDownMap = {}

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
        self.atoms[begin].nbr2.append(end)
        self.atoms[end].nbr2.append(begin)
        return bond

    def GetBond(self, idx):
        if idx >= 0 and idx < len(self.bonds):
            return self.bonds[idx]
        else:
            return None

    def NumBonds(self):
        return len(self.bonds)

    def GetBond(self, begin, end):
        for idx, bond in enumerate(self.bonds):
            if (bond.begin == begin and bond.end == end) or \
               (bond.begin == end and bond.end == begin):
                   return idx
        return -1

class Atom:
    def __init__(self):
        self.idx = -1
        self.element = ""
        self.charge = 0
        self.hcount = 0
        self.newidx = -1
        self.nbr = []
        self.nbr2 = []
        self.tetrahedralStereo = None
        self.ringClosure = False

    def __repr__(self):
        return "({}, idx: {}, newidx: {}, charge: {}, h: {}, nbr: {}, nbr2: {}, rc: {}, {})".format(self.element, self.idx, self.newidx, self.charge, self.hcount, self.nbr, self.nbr2, self.ringClosure, self.tetrahedralStereo)


class Bond:
    def __init__(self, begin=-1, end=-1, order=1):
        self.begin = begin
        self.end = end
        self.order = order

    def __repr__(self):
        return "(begin: {}, end: {}, order: {})".format(
                self.begin, self.end, self.order)


class ringClosure:
    def __init__(self, digit, prev, order, updown=''):
        self.digit = digit
        self.prev = prev
        self.order = order
        self.updown = updown

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
    updown = ''
    while stack:
        prev, idx = stack.pop()
        bondnum = mol.GetBond(prev, idx)
        if bondnum in mol.upDownMap:
            updown = mol.upDownMap[bondnum]
        if not visited[idx]:
            order.append(idx)
            mol.atoms[idx].newidx = len(order) - 1
            visited[idx] = True
            atom = mol.atoms[idx]
            # reorder tetrahedral stereochemistry
            if atom.tetrahedralStereo == 'Clockwise':
                atom.nbr[-2:] = [atom.nbr[-1], atom.nbr[-2]]
                atom.nbr2[-2:] = [atom.nbr2[-1], atom.nbr2[-2]]
            elif atom.tetrahedralStereo == 'CounterClockwise':
                pass

            # reorder cis/trans stereochemistry
            if mol.bonds[bondnum].order == 2:
                if len(atom.nbr) == 3:
                    bondnum = mol.GetBond(idx, atom.nbr[1])
                    bondnum2 = mol.GetBond(idx, atom.nbr[2])
                    if bondnum in mol.upDownMap and bondnum2 not in mol.upDownMap:
                        if mol.upDownMap[bondnum] == '\\':
                            mol.upDownMap[bondnum2] = '/'
                        elif mol.upDownMap[bondnum] == '/':
                            mol.upDownMap[bondnum2] = '\\'
                    if bondnum2 in mol.upDownMap and bondnum not in mol.upDownMap:
                        if mol.upDownMap[bondnum2] == '\\':
                            mol.upDownMap[bondnum] = '/'
                        elif mol.upDownMap[bondnum2] == '/':
                            mol.upDownMap[bondnum] = '\\'
                    if bondnum in mol.upDownMap and mol.upDownMap[bondnum] == updown:
                        atom.nbr[1], atom.nbr[2] = atom.nbr[2], atom.nbr[1]

            nbr = reversed(atom.nbr)
            for a in nbr:
                if not visited[a]:
                    stack.append((idx, a))
                elif a != prev:
                    mol.atoms[a].ringClosure = True
                    rcbond.append((a, idx))
    print(rcbond)

    for atom in mol.atoms:
        print(atom)
    #for bond in mol.bonds:
    #    print(bond)
    #print(mol.upDownMap)

    usedDigit = []
    rcstack = []
    pstack = []

    def generateBranchSmiles(idx):
        nonlocal mol, order, usedDigit, rcstack, rcbond, pstack
        smiles = ""
        updown = ''
        prev = pstack.pop()
        while True:
            atom = mol.atoms[idx]
            if atom.tetrahedralStereo is None:
                if atom.charge == 0 and atom.hcount == 0:
                    if atom.element == 'H':
                        smiles += '[H]'
                    else:
                        smiles += atom.element
                else:
                    s_charge = ''
                    if atom.charge > 0:
                        s_charge = '+'
                    elif atom.charge < 0:
                        s_charge = '-'
                    if abs(atom.charge) > 1:
                        s_charge += str(abs(atom.charge))
                    s_hcount = ''
                    if atom.hcount > 0:
                        s_hcount = 'H'
                        if atom.hcount > 1:
                            s_hcount += str(atom.hcount)
                    smiles += '[{}{}{}]'.format(atom.element, s_hcount, s_charge)
            else:
                smiles += '[' + atom.element + '@'
                if atom.hcount == 1:
                    smiles += 'H'
                smiles += ']'
            # if there is branch
            branch = [n for n in atom.nbr2 if n != prev]
            if not branch:
                return idx, smiles
            while len(branch) > 1:
                bidx = branch.pop(0)
                print(idx, bidx)
                # if ring is starting
                if (idx, bidx) in rcbond:
                    if not usedDigit:
                        digit = 1
                    else:
                        digit = max(usedDigit) + 1
                    usedDigit.append(digit)
                    rc = ringClosure(digit, idx, -1)
                    rcstack.append(rc)
                    smiles += str(digit)
                    continue
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
                    for i, b in enumerate(mol.bonds):
                        if (b.begin == idx and b.end == bidx) or (b.begin == bidx and b.end == idx):
                            if b.order == 1:
                                if i in mol.upDownMap:
                                    smiles += mol.upDownMap[i]
                            elif b.order == 2:
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
                isTransBranch = False
                for i, b in enumerate(mol.bonds):
                    if (b.begin == idx and b.end == nidx) or (b.begin == nidx and b.end == idx):
                        if b.order == 1:
                            if i in mol.upDownMap:
                                if updown == mol.upDownMap[i]:
                                    smiles += '('
                                    smiles += mol.upDownMap[i]
                                    pstack.append(idx)
                                    _, branch_smiles = generateBranchSmiles(nidx)
                                    smiles += branch_smiles
                                    smiles += ')'
                                    isTransBranch = True
                                    break
                                else:
                                    smiles += mol.upDownMap[i]
                                    updown = mol.upDownMap[i]
                        elif b.order == 2:
                            smiles += '='
                        elif b.order == 3:
                            smiles += '#'
                if isTransBranch:
                    return idx, smiles
                prev = idx
                idx = nidx
    pstack = [-1]
    _, smiles = generateBranchSmiles(0)
    return smiles

def validate(smiles):
    mol = smilesToMol(smiles)
    new_smiles = generateSmiles(mol)
    print(new_smiles)
    isCorrect = (Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) == Chem.MolToSmiles(Chem.MolFromSmiles(new_smiles)))
    return isCorrect, new_smiles

#with open('250k_rndm_zinc_drugs_clean.smi') as f:
#    for i, line in enumerate(f):
#        smiles = line.rstrip()
#        result, new_smiles = validate(smiles)
#        if not result:
#            print(i)
#            print(smiles)
#            print(new_smiles)
#            break

print(validate('CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C'))
