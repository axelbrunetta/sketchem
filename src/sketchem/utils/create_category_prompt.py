# This file contains the prompt for the Gemini AI to generate a new molecule category


def create_prompt(user_prompt):
    return f"""
Generate a list of molecule names that fit most accurately a category described by : "{user_prompt}".

Please provide 5-15 molecule names (except if a number was provided in the "text" from before, in which case use that one for the number of molecules).

Output Format:
Provide your response as a simple list of molecule names, each on a new line. Do NOT include SMILES notation, category names, the number of molecules, or any other text.

Category Name
Molecule 1 Name
Molecule 2 Name
...

For example, if the category was "Common Alcohols":
Common Alcohols
Methanol
Ethanol
Isopropanol

Or if the category was "Halogenated Solvents":
Halogenated Solvents
Chloroform
Dichloromethane

⸻

IMPORTANT Rules to Follow for Molecule Name Selection and Output:

1.  Names and Category Only: Your entire output should be the category name (without number of molecules) followed by a list of molecule names, each on a new line. Do not include ANY other explanations, commentary, or SMILES strings.
2.  Real Molecules Only: DO NOT invent molecules or use molecules that do not exist in real life (e.g., do NOT use ones that exist in movies or stories, or are fictional / imaginary). Use common, well-established chemical names.
3.  Do not add anything to the chemical names (e.g. don't write the name then something in parentheses next to it to justify the choice of that molecule)
4.  Simple Structures: Focus on molecules that can be easily represented by simple chemical structures (suitable for hand drawing). Avoid names that refer to:
    *   DNA, RNA
    *   Graphite, Diamond (extended covalent networks), Graphite, Graphene, ...
    *   Alloys
    *   Fullerenes, Carbon nanotubes, Buckyballs
    *   Complex / Very Large Proteins or Polymers
5.  Molecules, Not Atoms: Do NOT output names of elemental atoms only. Do NOT output things like:
    Hydrogen
    Helium
    Lithium
    Beryllium
    Boron
    Carbon
    Nitrogen
    Oxygen
    Fluorine
    Neon
    Sodium
    Magnesium
    Aluminum
    Silicon
    Phosphorus
    Sulfur
    Chlorine
    Argon
    Potassium
    Calcium
    Scandium
    Titanium
    Vanadium
    Chromium
    Manganese
    Iron
    Cobalt
    Nickel
    Copper
    Zinc
    Gallium
    Germanium
    Arsenic
    Selenium
    Bromine
    Krypton
    Rubidium
    Strontium
    Yttrium
    Zirconium
    Niobium
    Molybdenum
    Technetium
    Ruthenium
    Rhodium
    Palladium
    Silver
    Cadmium
    Indium
    Tin
    Antimony
    Tellurium
    Iodine
    Xenon
    Cesium
    Barium
    Lanthanum
    Cerium
    Praseodymium
    Neodymium
    Promethium
    Samarium
    Europium
    Gadolinium
    Terbium
    Dysprosium
    Holmium
    Erbium
    Thulium
    Ytterbium
    Lutetium
    Hafnium
    Tantalum
    Tungsten
    Rhenium
    Osmium
    Iridium
    Platinum
    Gold
    Mercury
    Thallium
    Lead
    Bismuth
    Polonium
    Astatine
    Radon
    Francium
    Radium
    Actinium
    Thorium
    Protactinium
    Uranium
    Neptunium
    Plutonium
    Americium
    Curium
    Berkelium
    Californium
    Einsteinium
    Fermium
    Mendelevium
    Nobelium
    Lawrencium
    Rutherfordium
    Dubnium
    Seaborgium
    Bohrium
    Hassium
    Meitnerium
    Darmstadtium
    Roentgenium
    Copernicium
    Nihonium
    Flerovium
    Moscovium
    Livermorium
    Tennessine
    Oganesson
6.  Category Interpretation: Be lenient on the category descriptions. If the description is vague, try to find molecule names related to that description (even if distantly related) while still respecting all aforementioned rules.
7.  Clarity of Names: Try to use names that are likely to be found in chemical databases like PubChem (e.g., common names or IUPAC names). Avoid overly obscure or ambiguous names.
"""