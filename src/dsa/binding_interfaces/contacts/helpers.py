import gemmi

#in PDB files some forms of lysine have special abbreviations
ALL_LYSINE_NAMES = {'LYS', 'LYZ', 'KCX', 'MLY', 'MLZ', 'M3L', 'ALY', 'ALZ'}

# def is_lysine(residue: gemmi.Residue|str) -> bool:
#     if isinstance(residue, gemmi.Residue):
#         return residue.name in ALL_LYSINE_NAMES
#     else:
#         return residue in ALL_LYSINE_NAMES

def is_lysine(residue: gemmi.Residue) -> bool:
    res_info = gemmi.find_tabulated_residue(residue.name)
    return res_info.one_letter_code.upper() == 'K'

def is_nz_atom(atom: gemmi.Atom|str) -> bool:
    if isinstance(atom, gemmi.Atom):
        return atom.name == 'NZ'
    else:
        return atom == 'NZ'

def has_nz_atom(residue: gemmi.Residue) -> bool:
    return any(is_nz_atom(atom) for atom in residue.atoms)


def is_lysine_and_nz(residue: gemmi.Residue, atom: gemmi.Atom) -> bool:
    return is_lysine(residue) and is_nz_atom(atom)

def is_atom_representative(residue: gemmi.Residue, atom: gemmi.Atom, polymer_type: gemmi.PolymerType) -> bool:
    if residue.entity_type != gemmi.EntityType.Polymer:
        return False
    if is_peptide(polymer_type):
        if residue.name == 'GLY':
            return atom.name == 'CA'
        else:
            return atom.name == 'CB'
    if is_dna_or_rna(polymer_type):
        return atom.name == "C1'"
    return False

def is_peptide(polymer_type: gemmi.PolymerType) -> bool:
    return polymer_type == gemmi.PolymerType.PeptideL or polymer_type == gemmi.PolymerType.PeptideD
def is_dna_or_rna(polymer_type: gemmi.PolymerType) -> bool:
    return polymer_type == gemmi.PolymerType.Dna or polymer_type == gemmi.PolymerType.Rna or polymer_type == gemmi.PolymerType.DnaRnaHybrid



def are_residues_adjacent(
    res1: gemmi.Residue,
    res2: gemmi.Residue,
    max_peptide_bond_dist: float = 1.5,
    max_nucleic_bond_dist: float = 1.8,
) -> bool:
    """Checks if res1 and res2 are sequentially and physically adjacent.

    Handles both proteins (C-N bond) and nucleic acids (O3'-P bond).
    """
    # 1. Verify sequence continuity (res2 must immediately follow res1)
    if (
        res1.label_seq is None
        or res2.label_seq is None
        or (res2.label_seq - res1.label_seq != 1)
    ):
        return False

    # 2. Check Physical Covalent Bond Distance
    res1_info = gemmi.find_tabulated_residue(res1.name)
    res2_info = gemmi.find_tabulated_residue(res2.name)
    # Case A: Polypeptides (Peptide Bond: C of res1 -> N of res2)
    if res1_info.is_amino_acid() and res2_info.is_amino_acid():
        c_atom = res1.find_atom("C", "\0")
        n_atom = res2.find_atom("N", "\0")

        if c_atom and n_atom:
            # dist_sq avoids computing square roots for speed
            return c_atom.pos.dist(n_atom.pos) < max_peptide_bond_dist

    # Case B: Polynucleotides (Phosphodiester Bond: O3' of res1 -> P of res2)
    elif res1_info.is_nucleic_acid() and res2_info.is_nucleic_acid():
        # Atom names for O3' can sometimes be written as O3* in older files
        o3_atom = res1.find_atom("O3'", "\0") or res1.find_atom("O3*", "\0")
        p_atom = res2.find_atom("P", "\0")

        if o3_atom and p_atom:
            return o3_atom.pos.dist(p_atom.pos) < max_nucleic_bond_dist

    return False

def characterize_residue(residue: gemmi.Residue) -> dict:
    """Characterizes a Gemmi residue into key biological and structural categories."""
    # Lookup built-in chemical component data for this residue name
    res_info = gemmi.find_tabulated_residue(residue.name)
    properties = {
        "is_amino_acid": res_info.is_amino_acid(),
        "is_nucleotide": res_info.is_nucleic_acid(),
        "is_water": res_info.is_water(),
        "is_standard": res_info.is_standard(),
        "is_polymer": residue.entity_type == gemmi.EntityType.Polymer,
        "is_lysine": res_info.one_letter_code.upper() == 'K',
    }
    is_ptm = properties['is_polymer'] and not properties['is_standard']
    is_ligand = (
        residue.entity_type == gemmi.EntityType.NonPolymer
        if residue.entity_type != gemmi.EntityType.Unknown
        else not (properties['is_polymer'] or properties['is_water'])
    )
    properties['is_ptm'] = is_ptm
    properties['is_ligand'] = is_ligand
    # has_ca = residue.get_ca() is not None
    # has_p = residue.get_p() is not None
    # is_polymer = (
    #     residue.entity_type == gemmi.EntityType.Polymer or has_ca or has_p
    # )
    # is_amino_acid = res_info.is_amino_acid() or has_ca
    # is_nucleotide = res_info.is_nucleic_acid() or (has_p and not has_ca)
    # is_water = res_info.is_water()
    # is_standard = res_info.is_standard()
    # is_ptm = is_polymer and not is_standard
    # name_upper = residue.name.upper()
    # is_lysine = name_upper in ("LYS", "K") or (
    #     is_ptm and name_upper in ("MLY", "KCX", "MLZ", "M3L")
    # )
    return properties


