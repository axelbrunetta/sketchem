MOLECULE_CATEGORIES = {
    "Alkanes (8)": {
        "Methane": "C",
        "Ethane": "CC",
        "Propane": "CCC",
        "Butane": "CCCC",
        "Isobutane": "CC(C)C",
        "Pentane": "CCCCC",
        "Isopentane": "CCC(C)C",
        "Neopentane": "CC(C)(C)C",
    },
    "Alkenes (6)": {
        "Ethene": "C=C",
        "Propene": "CC=C",
        "1-Butene": "CCC=C",
        "2-Butene": "C/C=C/C",
        "Isobutene": "CC(=C)C",
        "1,3-Butadiene": "C=CC=C",
    },
    "Alkynes (6)": {
        "Ethyne (Acetylene)": "C#C",
        "Propyne": "CC#C",
        "1-Butyne": "CCC#C",
        "2-Butyne": "CC#CC",
        "1-Pentyne": "CCCC#C",
        "3-Methyl-1-butyne": "CC(C)C#C",
    },
    "Common Solvents (15)": {
        "Water": "O",
        "Ethanol": "CCO",
        "Methanol": "CO",
        "Isopropanol": "CC(O)C",
        "Acetone": "CC(=O)C",
        "Ethyl acetate": "CC(=O)OCC",
        "Diethyl ether": "CCOCC",
        "Tetrahydrofuran (THF)": "C1CCOC1",
        "Dichloromethane (DCM)": "ClCCl",
        "Chloroform": "ClC(Cl)Cl",
        "Hexane": "CCCCCC",
        "Heptane": "CCCCCCC",
        "Toluene": "Cc1ccccc1",
        "Dimethylformamide (DMF)": "CN(C)C=O",
        "Dimethyl sulfoxide (DMSO)": "CS(=O)C",
    },
    "Carboxylic Acids (8)": {
        "Formic acid": "C(=O)O",
        "Acetic acid": "CC(=O)O",
        "Propanoic acid": "CCC(=O)O",
        "Butanoic acid": "CCCC(=O)O",
        "Benzoic acid": "C(=O)(O)c1ccccc1",
        "Oxalic acid": "C(=O)(C(=O)O)O",
        "Malonic acid": "C(C(=O)O)C(=O)O",
        "Succinic acid": "C(CC(=O)O)C(=O)O",
    },
    
    "Polycyclic Aromatic Hydrocarbons (PAHs) (6)": {
        "Naphthalene": "c1ccc2ccccc2c1",
        "Anthracene": "c1ccc2cc3ccccc3cc2c1",
        "Phenanthrene": "c1ccc2c(c1)ccc3ccccc23",
        "Chrysene": "c1ccc2c(c1)ccc3c2ccc4ccccc34",
        "Coronene": "c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67",
        "Benzo[a]pyrene": "c1ccc2c(c1)cc3ccc4cccc5c4c3c2C=C5"
    },
    "Aromatic Heterocycles (8)": {
        "Quinoline": "c1ccc2ncccc2c1",
        "Isoquinoline": "c1ccc2cnccc2c1",
        "Indole": "C12=C(C=CN2)C=CC=C1",
        "Benzofuran": "c1ccc2c(c1)occ2",
        "Benzothiophene": "c1ccc2c(c1)scc2",
        "Carbazole": "c1ccc2c(c1)c3ccccc3[nH]2",
        "Acridine": "c1ccc2nc3ccccc3cc2c1",
        "Purine": "c1c2c(nc[nH]2)ncn1"
    },
    "Amino Acids (L-isomers) (9)": {
        "L-Leucine": "CC(C)C[C@H](N)C(=O)O",
        "L-Isoleucine": "CC[C@H](C)[C@H](N)C(=O)O",
        "L-Valine": "CC(C)[C@H](N)C(=O)O",
        "L-Phenylalanine": "N[C@H](C(=O)O)Cc1ccccc1",
        "L-Tryptophan": "N[C@H](C(=O)O)Cc1c[nH]c2ccccc12",
        "L-Methionine": "CSCC[C@H](N)C(=O)O",
        "L-Threonine": "C[C@@H](O)[C@H](N)C(O)=O",
        "L-Lysine": "NCCCC[C@H](N)C(=O)O",
        "L-Histidine": "N[C@H](C(=O)O)Cc1c[nH]cn1"
    },
    "Common Monosaccharides (Cyclic,D-isomers) (6)": {
        "Alpha-D-Glucopyranose": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Beta-D-Glucopyranose": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "Alpha-D-Galactopyranose": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
        "Alpha-D-Mannopyranose": "OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
        "Beta-D-Fructofuranose": "C([C@@H]1[C@H]([C@@H]([C@](O1)(CO)O)O)O)O",
        "Alpha-D-Ribofuranose": "C([C@@H]1[C@H]([C@H]([C@H](O1)O)O)O)O"
    },
    "Nucleobases & Common Nucleosides (6)": {
        "Adenine": "Nc1ncnc2[nH]cnc12",
        "Guanine": "C1=NC2=C(N1)C(=O)NC(=N2)N",
        "Cytosine": "C1=C(NC(=O)N=C1)N", 
        "Thymine": "CC1=CNC(=O)NC1=O", 
        "Uracil": "O=c1ccn([H])c(=O)[nH]1", 
        "Adenosine": "Nc1ncnc2n(cnc12)[C@@H]3O[C@H](CO)[C@@H](O)[C@H]3O"
    },
    "Common Small Pharmaceuticals (7)": {
        "Aspirin (Acetylsalicylic acid)": "CC(=O)Oc1ccccc1C(=O)O",
        "Paracetamol (Acetaminophen)": "CC(=O)Nc1ccc(O)cc1",
        "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "Metformin": "CN(C)C(=N)N=C(N)N",
        "Diazepam (Valium)": "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3"
    },

}

