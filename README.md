# Stereo SMILES Converter
```
$ python parser.py 
True F/C=C(\Cl)/I -> F/C=C(\Cl)/I
True F/C=C(/Cl)\I -> F/C=C(\I)/Cl
True F/C=C/I -> F/C=C([H])/I
True F/C=C\I -> F/C=C(\I)[H]
True F/C=C/F -> F/C=C([H])/F
True N[C@](O)(Br)C -> N[C@](O)(Br)C
True N[C@@](O)(Br)C -> N[C@](O)(C)Br
True N[C@H](O)C -> N[C@H](O)C
True N[C@@H](O)C -> N[C@H](C)O
True [C@H](N)(O)C -> [C@H](N)(O)C
True [C@@H](N)(O)C -> [C@H](N)(C)O
True [C@](Br)(N)(O)C -> [C@](Br)(N)(O)C
True [C@@](Br)(N)(O)C -> [C@](Br)(N)(C)O
True [C@]1(Br)(Cl)CCCC(F)C1 -> [C@]1(Br)(Cl)CCCC(F)C1
True [C@@]1(Br)(Cl)CCCC(F)C1 -> [C@]1(Br)(Cl)CC(CCC1)F
True N1C[C@H]2CC=CC[C@@H]2C1 -> N1C[C@H]2CC=CC[C@H](C1)2
True C[C@@H]1CC(Nc2cncc(-c3nncn3C)c2)C[C@@H](C)C1 -> C[C@H]1C[C@H](C)CC(C1)Nc1cncc(c2nncn2C)c1
True N#Cc1ccc(-c2ccc(O[C@@H](C(=O)N3CCCC3)c3ccccc3)cc2)cc1 -> N#Cc1ccc(c2ccc(O[C@H](c3ccccc3)C(=O)N3CCCC3)cc2)cc1
True CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C1 -> CCOC(=O)[C@H]1CN(CCC1)C(=O)c1nc(c2ccc(C)cc2)n2c1CCCCC2
True N#CC1=C(SCC(=O)Nc2cccc(Cl)c2)N=C([O-])[C@H](C#N)C12CCCCC2 -> N#CC1=C(SCC(=O)Nc2cccc(Cl)c2)N=C([O-])[C@H](C#N)C21CCCCC2
True CC[NH+](CC)[C@](C)(CC)[C@H](O)c1cscc1Br -> CC[NH+](CC)[C@](C)(CC)[C@H](O)c1cscc1Br
True C[C@@H]1CN(C(=O)c2cc(Br)cn2C)CC[C@H]1[NH3+] -> C[C@H]1[C@H]([NH3+])CCN(C1)C(=O)c1cc(Br)cn1C
True CCOc1ccc(OCC)c([C@H]2C(C#N)=C(N)N(c3ccccc3C(F)(F)F)C3=C2C(=O)CCC3)c1 -> CCOc1ccc(OCC)c([C@H]2C(C#N)=C(N)N(c3ccccc3C(F)(F)F)C4=C2C(=O)CCC4)c1
True Cc1ccc2nc(S[C@H](C)C(=O)NC3CCC(C)CC3)n(C)c(=O)c2c1 -> Cc1ccc2nc(S[C@H](C)C(=O)NC3CCC(C)CC3)n(C)c(=O)c2c1
True Cc1ccccc1C(=O)N1CCC2(CC1)C[C@H](c1ccccc1)C(=O)N2C -> Cc1ccccc1C(=O)N2CCC3(CC2)C[C@H](c4ccccc4)C(=O)N3C
True CC(C)Cc1nc(SCC(=O)NC[C@@H]2CCCO2)c2c(=O)n(C)c(=O)n(C)c2n1 -> CC(C)Cc1nc(SCC(=O)NC[C@H]2OCCC2)c2c(=O)n(C)c(=O)n(C)c2n1
True Cc1ccc(CNC(=O)c2ccccc2NC(=O)[C@@H]2CC(=O)N(c3ccc(C)cc3)C2)cc1 -> Cc1ccc(CNC(=O)c2ccccc2NC(=O)[C@H]3CN(C(C3)=O)c3ccc(C)cc3)cc1
True CC(C)[C@@H](Oc1cccc(Cl)c1)C(=O)N1CCC(n2cccn2)CC1 -> CC(C)[C@H](C(=O)N1CCC(n2cccn2)CC1)Oc1cccc(Cl)c1
True CCN(CC)C(=O)C[C@@H](C)[NH2+][C@H](C)c1cccc(F)c1F -> CCN(CC)C(=O)C[C@H]([NH2+][C@H](C)c1cccc(F)c1F)C
True O=C(NCCNC(=O)N1C[C@H]2CC=CC[C@@H]2C1)c1cccnc1 -> O=C(NCCNC(=O)N1C[C@H]2CC=CC[C@H](C1)2)c1cccnc1
True Cc1ccc(N2CC[C@@H](NS(=O)(=O)c3ccccc3C)C2=O)cc1C -> Cc1ccc(N2CC[C@H](C2=O)NS(=O)(=O)c3ccccc3C)cc1C
True CC[C@H](C)C[C@@H](C)NC(=O)N1CCN(CC(=O)NC2CC2)CC1 -> CC[C@H](C)C[C@H](NC(=O)N1CCN(CC(=O)NC2CC2)CC1)C
True CC(=O)Nc1c2n(c3ccccc13)C[C@](C)(C(=O)NC1CCCCC1)N(C1CCCCC1)C2=O -> CC(=O)Nc1c2n(c3ccccc13)C[C@](C)(C(=O)NC3CCCCC3)N(C3CCCCC3)C2=O
True N#Cc1ccncc1NC[C@@H]1C[C@@]12CCc1ccccc12 -> N#Cc1ccncc1NC[C@H]2[C@]3(C2)c4ccccc(CC3)4
True CCOC(=O)c1nnc2ccccc2c1N1CC[C@@H]([NH+](CC)CC)C1 -> CCOC(=O)c1nnc2ccccc2c1N3CC[C@H](C3)[NH+](CC)CC
True CCC[NH2+][C@@H]1COC[C@H]1C(=O)NCc1cscc1C -> CCC[NH2+][C@H]1[C@H](C(=O)NCc2cscc2C)COC1
True CC(=O)c1ccc(S(=O)(=O)N2CCCC[C@H]2C)cc1 -> CC(=O)c1ccc(S(=O)(=O)N2CCCC[C@H]2C)cc1
```
