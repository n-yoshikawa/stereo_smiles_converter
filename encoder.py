class smilesEncoder:
    def __init__(self):
        self.pos = 0

    def encode_smiles(self, smiles):
        self.pos = 0
        return self.get_branch(smiles)


    def get_branch(self, smiles):
        chain1 = self.get_chain(smiles)
        if self.pos == len(smiles) or smiles[self.pos] == ')':
            return chain1
        self.pos += 1  # skip '['
        bracket_atom = ''
        while smiles[self.pos] != ']':
            bracket_atom += smiles[self.pos]
            self.pos += 1
        clockwise = False
        if '@@' in bracket_atom:
            clockwise = True
        bracket_atom = bracket_atom.replace('@', '')
        bracket_atom = bracket_atom.replace('H', '')
        self.pos += 1  # skip ']'
        branch1 = ''
        if smiles[self.pos] == '(':
            self.pos += 1  # skip '('
            branch1 = self.get_branch(smiles)
            self.pos += 1  # skip ')'
        branch2 = ""
        if smiles[self.pos] == '(':
            self.pos += 1  # skip '('
            branch2 = self.get_branch(smiles)
            self.pos += 1  # skip ')'
        chain2 = self.get_branch(smiles)
        if clockwise:
            if branch2 == '':
                return "{}{}({}){}".format(chain1, bracket_atom, chain2, branch1)
            else:
                return "{}{}({})({}){}".format(chain1, bracket_atom,
                                               branch2, branch1, chain2)
        else:
            if branch2 == '':
                return "{}{}({}){}".format(chain1, bracket_atom, branch1, chain2)
            else:
                return "{}{}({})({}){}".format(chain1, bracket_atom,
                                               branch1, branch2, chain2)

    def get_chain(self, smiles):
        chain = ''
        while self.pos < len(smiles) and smiles[self.pos] not in ['[', ')']:
            chain += smiles[self.pos]
            self.pos += 1
        return chain


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
print("Unsupported case:")
smiles = "[C@@](Br)(N)(O)C"
print(smiles, '->', a.encode_smiles(smiles))
smiles = "[C@]1(Br)(Cl)CCCC(F)C1"
print(smiles, '->', a.encode_smiles(smiles))
