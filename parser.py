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
        if element.isupper() and ptr+1 < len(smiles) and smiles[ptr+1].islower():
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
    usedDigit = []
    rcstack = []

    # decide atom order and find ring closure
    order = []
    visited = [False for _ in range(len(mol.atoms))]
    stack = [(-1, 0)]
    while stack:
        prev, idx = stack.pop()
        if visited[idx] and prev > idx+1:
            mol.atoms[idx].ringClosure = True
        if not visited[idx]:
            order.append(idx)
            mol.atoms[idx].newidx = len(order) - 1
            visited[idx] = True
            atom = mol.atoms[idx]
            nbr = atom.nbr
            if atom.tetrahedralStereo != 'Clockwise':
                nbr = reversed(nbr)
            for a in nbr:
                stack.append((idx, a))

    def generateBranchSmiles(idx):
        nonlocal mol, order, usedDigit, rcstack
        smiles = ""
        while True:
            atom = mol.atoms[idx]
            if atom.tetrahedralStereo is None:
                smiles += atom.element
            else:
                smiles += '[' + atom.element + '@H]'
            # if ring is starting
            if atom.ringClosure:
                if not usedDigit:
                    digit = 1
                else:
                    digit = min(usedDigit) + 1
                usedDigit.append(digit)
                rc = ringClosure(digit, idx, -1)
                rcstack.append(rc)
                smiles += str(digit)
            # if branching (buggy)
            numBranch = len([mol.atoms[n].newidx > atom.newidx for n in atom.nbr])
            if (atom.ringClosure and numBranch > 1) or (not atom.ringClosure and numBranch > 2):
                smiles += '('
                idx, branch = generateBranchSmiles(idx+1)
                smiles += branch
                smiles += ')'
            ringDest = [mol.atoms[n].newidx for n in atom.nbr if mol.atoms[n].newidx < atom.newidx]
            for rc in rcstack:
                if rc.prev < atom.newidx - 1 and rc.prev in ringDest:
                    smiles += str(rc.digit)
                    usedDigit.remove(rc.digit)
            if idx+1 >= len(order) or order[idx+1] not in atom.nbr:
                return idx, smiles
            else:
                for b in mol.bonds:
                    if b.begin == idx and b.end == idx + 1:
                        if b.order == 2:
                            smiles += '='
                        elif b.order == 3:
                            smiles += '#'
                idx += 1
    _, smiles = generateBranchSmiles(0)
    return smiles

def evaluate(smiles):
    mol = smilesToMol(smiles)
    print(smiles, '->', generateSmiles(mol))

evaluate('N1C[C@H]2CC=CC[C@@H]2C1')
