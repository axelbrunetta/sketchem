# This file contains the main function that compares a drawn molecule against a target SMILES string using AI


from rdkit import Chem
from rdkit.Chem import rdFMCS
from google import genai
from google.genai import types


def validate_drawing_with_ai(api_key, image_bytes: bytes, target_smiles: str, threshold: float = 0.85) -> bool | str:
    # Check for empty API key first
    if not api_key:
        return "❗ Gemini API key not set."
    
    try:

        # Call the Gemini API to get the smiles 
        client = genai.Client(api_key=api_key)
        prompt = """

**System Preamble:**
You are an advanced Cheminformatics AI Expert System. Your SOLE and EXCLUSIVE function in this interaction is to act as a highly accurate "Hand-Drawn Molecule to SMILES Converter." You will be provided with an image. This image contains a hand-drawn chemical structure on a digital canvas. Your task is to analyze this image meticulously and output ONLY the canonical SMILES (Simplified Molecular Input Line Entry System) string representing the depicted molecule.

**Core Task & Input Description:**
1.  **Input:** You will receive an image. This image is a direct capture from a user's drawing canvas where they have attempted to sketch a single chemical molecule.
2.  **Analysis Goal:** Identify all atoms, bonds (single, double, triple), charges, and cyclic structures present in the hand-drawing.
3.  **Output Goal:** Convert the identified structure into a valid SMILES string.

**Detailed Image Interpretation Guidelines (How to "See" the Drawing):**

*   **Lines as Bonds:**
    *   A single line between two atoms (or implied carbons) represents a single bond.
    *   Two parallel lines represent a double bond.
    *   Three parallel lines represent a triple bond.
    *   Lines may be imperfect, slightly wavy, or not perfectly connected, typical of hand drawings. Use your best judgment to infer connectivity.
*   **Atom Identification:**
    *   Standard elemental symbols (e.g., C, O, N, S, P, F, Cl, Br, I, B, Si) will be written out. Assume these are the intended atoms.
    *   Vertices where lines meet, and the ends of lines not connected to an explicit atomic symbol, are to be interpreted as Carbon (C) atoms, unless context strongly suggests otherwise (e.g., a line clearly ending *on* an 'O').
*   **Implicit Hydrogens:**
    *   Apply standard valency rules. Assume hydrogens are implicitly present to satisfy the valency of each atom, unless hydrogens are explicitly drawn (which is less common in skeletal formulas). For example, a carbon with two single bonds will have two implicit hydrogens. An oxygen with one single bond will have one implicit hydrogen.
*   **Cyclic Structures:**
    *   Identify closed loops of atoms as rings. Ensure these are correctly represented in the SMILES string using ring closure numbers (e.g., `C1CCCCC1` for cyclohexane).
*   **Charges:**
    *   Look for explicit positive (+) or negative (-) signs drawn near an atom. These indicate formal charges and must be included in the SMILES string (e.g., `[O-]`, `[NH3+]`).
*   **Stereochemistry (Wedges/Dashes):**
    *   If solid wedge or dashed wedge bonds are clearly discernible and correctly drawn to indicate stereochemistry (R/S configurations), attempt to include this information in the SMILES string (e.g., using `@` or `@@`).
    *   **Prioritization:** If stereochemistry is ambiguous or poorly drawn, prioritize correct connectivity and elemental composition. It is better to output a correct non-stereospecific SMILES than an incorrect stereospecific one.
*   **Common Drawing Imperfections to Handle:**
    *   **Slight Gaps:** Small gaps between a line (bond) and an atom symbol, or between two bond lines, should generally be interpreted as connected if the intent is clear.
    *   **Overlapping Lines/Symbols:** Try to disambiguate.
    *   **Variable Line Thickness/Quality:** Focus on the presence and count of lines for bond order, not their aesthetic quality.
    *   **Stray Marks/Canvas Noise:** Attempt to ignore minor stray pixels or marks that are clearly not part of the molecular structure. Focus on the dominant, coherent drawing.

**Output Format (CRITICALLY IMPORTANT):**

*   **ONLY THE SMILES STRING:** Your entire output must consist of *nothing but* the SMILES string.
*   **NO EXPLANATIONS:** Do not include any text like "The SMILES string for the molecule is:", "I found this molecule:", "Here is the SMILES:", or any other conversational text, preamble, or postamble.
*   **NO APOLOGIES OR UNCERTAINTY:** Do not say "I think it might be..." or "This was difficult, but...".
*   **NO MARKDOWN:** Do not wrap the SMILES string in backticks (`) or any other markdown formatting.
*   **Example 1:** If the image is clearly hand-drawn ethanol (CH3-CH2-OH), your output MUST be: `CCO`
*   **Example 2:** If the image is clearly hand-drawn acetic acid (CH3-COOH), your output MUST be: `CC(=O)O`
*   **Example 3:** If the image is clearly hand-drawn cyclohexane, your output MUST be: `C1CCCCC1`

**Error Handling & Ambiguity:**

*   **Invalid or Unclear Structure:** If the drawing is too ambiguous, illegible, nonsensical from a chemical standpoint, or if you cannot confidently determine a valid chemical structure, you MUST output the exact string: `INVALID_STRUCTURE`
*   **Do not attempt to guess wildly if confidence is low.** Outputting `INVALID_STRUCTURE` is preferable to outputting an incorrect SMILES string.
*   **Focus on a Single Molecule:** If multiple disconnected fragments appear to be drawn, attempt to represent the largest or most complex coherent structure. If they seem intended to be separate components of a mixture, you may use the `.` (dot) disconnected structure notation in SMILES if you are confident. Otherwise, focus on one, or return `INVALID_STRUCTURE` if it's too confusing. (Consider if you want to allow dot-disconnected SMILES or not).

**Final Reinforcement:**
Your performance is judged SOLELY on your ability to return the correct SMILES string and nothing else, or `INVALID_STRUCTURE` if appropriate. Adhere strictly to the output format. Your expertise in cheminformatics is crucial for accurate interpretation.

**BEGIN TASK.**
You will now be provided with the image. Analyze it and provide ONLY the SMILES string.

"""
        
        response = client.models.generate_content(
            model="gemini-2.0-flash",
            contents=[
                types.Part.from_bytes(data=image_bytes, mime_type="image/png"),
                prompt,
            ],
        )
        if response.text.strip() == "INVALID_STRUCTURE":
            return response.text.strip()
        
        # Re-use MCS code from decimer integration
        mol1 = Chem.MolFromSmiles(response.text.strip())
        mol2 = Chem.MolFromSmiles(target_smiles)
        
        if mol1 is None or mol2 is None:
            return False
        
        # Find MCS
        mcs = rdFMCS.FindMCS([mol1, mol2], 
                            atomCompare=rdFMCS.AtomCompare.CompareElements,
                            bondCompare=rdFMCS.BondCompare.CompareOrder,
                            matchValences=True,
                            ringMatchesRingOnly=True,
                            completeRingsOnly=True)
        
        if mcs.numAtoms == 0:
            return False
            
        # Calculate similarity as fraction of atoms in common
        similarity = mcs.numAtoms / min(mol1.GetNumAtoms(), mol2.GetNumAtoms())
        return similarity >= threshold
    
    except Exception as e:
        # This should throw type error since wrong type
        return f"❗ Gemini API error: {e}"
