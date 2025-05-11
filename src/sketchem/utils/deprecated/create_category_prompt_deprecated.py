def create_prompt(user_prompt):
        return f"""
Generate a list of molecules that fit most accurately a category described by : "{user_prompt}".

Please provide 5-10 molecules (except if a number was provided in the "text" from before, in which case use that one for the number of molecules).

**SMILES Notation Guidelines (CRITICALLY IMPORTANT):**
When providing the SMILES notation, you MUST adhere to the following to ensure correctness and avoid common errors:

1.  **SMILES for Connectivity, Not Just Formulas:**
    *   SMILES strings describe atomic connectivity. They are NOT simply molecular formulas.
    *   Numbers in SMILES (unless part of an element symbol like `[13C]`) are typically for ring closures, not atomic counts.

2.  **Common Pitfalls & Correct SMILES Examples:**
    *   **Methane (CH₄):**
        *   Correct SMILES: `C` (Hydrogens are implicit)
        *   Incorrect: `CH4` (The '4' would be misinterpreted as a ring number)
    *   **Water (H₂O):**
        *   Correct SMILES: `O` (Hydrogens are implicit)
        *   Incorrect: `H2O`
    *   **Ammonia (NH₃):**
        *   Correct SMILES: `N` (Hydrogens are implicit)
        *   Incorrect: `NH3`
    *   **Boron Trifluoride (BF₃):**
        *   Correct SMILES: `B(F)(F)F` or `FB(F)F` (Shows B connected to three F atoms)
        *   Incorrect: `BF3`
    *   **Borane (BH₃):**
        *   Correct SMILES: `B` (Hydrogens implicit) or `[BH3]` (Explicit hydrogens)
        *   Incorrect: `BH3`
    *   **Sulfur Hexafluoride (SF₆):**
        *   Correct SMILES: `FS(F)(F)(F)(F)F`
        *   Incorrect: `SF6`
3. ** More Correct SMILES Examples:**
*   **Simple Alkanes & Connectivity:**
    *   Ethane (CH3-CH3): `CC`
    *   Propane (CH3-CH2-CH3): `CCC`
    *   Isobutane (CH3-CH(CH3)-CH3): `CC(C)C`
    *   Neopentane (C(CH3)4): `CC(C)(C)C`

*   **Alkenes & Alkynes (Double & Triple Bonds):**
    *   Ethene (CH2=CH2): `C=C`
    *   Propene (CH3-CH=CH2): `CC=C`
    *   2-Butene (CH3-CH=CH-CH3): `CC=CC` (cis/trans ignored if not drawn)
    *   Ethyne (HC≡CH): `C#C`
    *   Propyne (CH3-C≡CH): `CC#C`
    *   1,3-Butadiene (CH2=CH-CH=CH2): `C=CC=C`

*   **Alcohols, Ethers, Thiols, Thioethers:**
    *   Methanol (CH3-OH): `CO`
    *   Ethanol (CH3-CH2-OH): `CCO`
    *   Isopropanol (CH3-CH(OH)-CH3): `CC(O)C`
    *   Dimethyl ether (CH3-O-CH3): `COC`
    *   Diethyl ether (CH3CH2-O-CH2CH3): `CCOCC`
    *   Ethanethiol (CH3-CH2-SH): `CCS`
    *   Dimethyl sulfide (CH3-S-CH3): `CSC`

*   **Aldehydes & Ketones:**
    *   Formaldehyde (CH2=O): `C=O`
    *   Acetaldehyde (CH3-CHO): `CC=O`
    *   Acetone (CH3-CO-CH3): `CC(=O)C`
    *   Butanone (CH3-CO-CH2-CH3): `CCC(=O)C`

*   **Carboxylic Acids & Esters:**
    *   Formic acid (HCOOH): `C(=O)O`
    *   Acetic acid (CH3-COOH): `CC(=O)O`
    *   Propionic acid (CH3-CH2-COOH): `CCC(=O)O`
    *   Methyl acetate (CH3-COO-CH3): `COC(=O)C` (Note order: `C-O-C=O` for ester linkage)
    *   Ethyl acetate (CH3-COO-CH2CH3): `CCOC(=O)C`

*   **Amines & Amides:**
    *   Methylamine (CH3-NH2): `CN`
    *   Dimethylamine ((CH3)2-NH): `CNC`
    *   Trimethylamine ((CH3)3-N): `CN(C)C`
    *   Ethylamine (CH3-CH2-NH2): `CCN`
    *   Acetamide (CH3-CO-NH2): `CC(=O)N`
    *   N-Methylacetamide (CH3-CO-NH-CH3): `CNC(=O)C`

*   **Halogenated Compounds:**
    *   Chloromethane (CH3-Cl): `CCl`
    *   Dichloromethane (CH2Cl2): `ClCCl` (or `C(Cl)Cl`)
    *   Chloroform (CHCl3): `ClC(Cl)Cl`
    *   Carbon tetrachloride (CCl4): `ClC(Cl)(Cl)Cl`
    *   Bromoethane (CH3-CH2-Br): `CCBr`
    *   Trifluoromethyl group (-CF3): `C(F)(F)F` (as part of a larger molecule, e.g., `CC(F)(F)F`)

*   **Molecular Formulas vs. SMILES (CRITICAL REINFORCEMENT):**
    *   Water (H2O, drawn as O with two H): `O` (NOT H2O)
    *   Ammonia (NH3, drawn as N with three H or just N): `N` (NOT NH3)
    *   Boron Trifluoride (BF3, drawn as B bonded to 3 F): `B(F)(F)F` (or `FB(F)F`) (NOT BF3)
    *   Borane (BH3, drawn as B with three H or just B): `B` (or `[BH3]`) (NOT BH3)
    *   Silicon Tetrafluoride (SiF4): `F[Si](F)(F)F` (NOT SiF4)
    *   Sulfur Hexafluoride (SF6): `FS(F)(F)(F)(F)F` (NOT SF6)

*   **Cyclic Structures (Aliphatic):**
    *   Cyclopropane: `C1CC1`
    *   Cyclobutane: `C1CCC1`
    *   Cyclopentane: `C1CCCC1`
    *   Cyclohexane: `C1CCCCC1`
    *   Methylcyclohexane: `CC1CCCCC1`
    *   1,1-Dimethylcyclopropane: `CC1(C)CC1`
    *   Cyclohexene: `C1CC=CCC1`

*   **Aromatic Systems (Benzene & Derivatives):**
    *   Benzene: `c1ccccc1` (or `C1=CC=CC=C1` - Kekule form, but aromatic form preferred if recognized)
    *   Toluene (Methylbenzene): `Cc1ccccc1`
    *   Phenol (Hydroxybenzene): `Oc1ccccc1`
    *   Aniline (Aminobenzene): `Nc1ccccc1`
    *   Nitrobenzene: `O=[N+]([O-])c1ccccc1` (or `[O-][N+](=O)c1ccccc1`)
    *   Chlorobenzene: `Clc1ccccc1`
    *   Styrene (Vinylbenzene): `C=Cc1ccccc1`
    *   Naphthalene (fused rings): `c1ccc2ccccc2c1`
    *   Anthracene: `c1ccc2cc3ccccc3cc2c1`

*   **Heterocyclic Aromatic Systems:**
    *   Pyridine: `n1ccccc1`
    *   Pyrrole: `c1cc[nH]c1`
    *   Furan: `c1ccoc1`
    *   Thiophene: `c1ccsc1`
    *   Imidazole: `c1c[nH]cn1`
    *   Quinoline: `n1ccc2ccccc2c1`

4.  **Branching and Bonds:**
    *   Use parentheses `()` for branches: e.g., Isobutane `CC(C)C`
    *   Use `=` for double bonds: e.g., Ethene `C=C`
    *   Use `#` for triple bonds: e.g., Ethyne `C#C`

5.  **Cyclic Structures:**
    *   Use ring closure numbers: e.g., Cyclohexane `C1CCCCC1`
    *   Aromatic rings often use lowercase: e.g., Benzene `c1ccccc1` (or `C1=CC=CC=C1` in Kekule form)

6.  **Heteroatoms and Functional Groups:**
    *   Ethanol (CH₃CH₂OH): `CCO`
    *   Acetic Acid (CH₃COOH): `CC(=O)O`
    *   Methylamine (CH₃NH₂): `CN`

**Output Format:**
Provide your response in the following format ONLY:
Category Name (number of molecules)
Molecule 1 Name: SMILES notation
Molecule 2 Name: SMILES notation
...

**For example (Note the correct SMILES usage):**
Common Alcohols (3)
Methanol: CO
Ethanol: CCO
Isopropanol: CC(O)C

Simple Inorganic Hydrides (2)
Water: O
Ammonia: N

Fluorinated Boron Compounds (1)
Boron Trifluoride: B(F)(F)F
⸻

**IMPORTANT Rules to Follow for Molecule Selection and Output:**

1.  **No Extra Commentary:** Before and after providing this formatting of the category name, molecule names, and their SMILES, do not include ANY other explanations or commentary. Simply output what is asked above.
2.  **Real Molecules Only:** DO NOT invent molecules or use molecules that do not exist in real life (e.g., do NOT use ones that exist in movies or stories, or are fictional / imaginary).
3.  **Simple Structures:** Make sure the molecules you find can be easily represented using a simple chemical structure drawing with bonds and atoms. Avoid complex structures like:
    *   DNA, RNA
    *   Graphite, Diamond (extended covalent networks)
    *   Alloys
    *   Fullerenes, Carbon nanotubes
    *   Complex / Very Large Proteins or Polymers
4.  **Molecules, Not Atoms:** Do NOT output elemental atoms only (I am looking for molecules). Do NOT output things like:
    *   Osmium
    *   Bismuth
    *   Radon
    *   Krypton
    *   Iron
    *   Gold
5.  **Valid SMILES:** Ensure the SMILES you generate are valid according to the guidelines above. If you are not sure about a molecule's correct SMILES representation or if it violates the guidelines, try to find another molecule that fits the category and has a clearly representable, valid SMILES string.
6.  **Category Interpretation:** Be lenient on the category descriptions. If the description is vague, try to find molecules related to that description (even if distantly related) while still respecting all aforementioned rules, especially the SMILES notation guidelines.

"""