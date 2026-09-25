import requests
import os
from dsa.uniprot.config import ALL_UNIPROT_IDS
from dsa.binding_interfaces.config import STRUCTURE_DIR


def download_alphafold_structure(uniprot_id, out_dir=".", version=4, format="pdb"):
    """
    Download AlphaFold structure for a given UniProt ID.
    
    Args:
        uniprot_id : UniProt accession (e.g. 'P12345')
        out_dir    : output directory
        version    : AlphaFold model version (default 4)
        format     : 'pdb' or 'cif'
    """
    if format != "pdb" and format != "cif":
        raise ValueError("format must be 'pdb' or 'cif'")
    file_name = f"AF-{uniprot_id}-F1-model_v{version}.{format}"
    # if format == "pdb":
    #     url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.pdb"
    # elif format == "cif":
    #     url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.cif"
    # else:
    #     raise ValueError("format must be 'pdb' or 'cif'")
    url = f"https://alphafold.ebi.ac.uk/files/{file_name}"

    out_path = out_dir / file_name

    response = requests.get(url)
    if response.status_code == 200:
        # with open(out_path, "wb") as f:
        #     f.write(response.content)
        # print(f"Downloaded: {out_path}")
        return out_path
    elif response.status_code == 404:
        print(f"Not found: {uniprot_id} (no AlphaFold model available)")
        return None
    else:
        print(f"Error {response.status_code} for {uniprot_id}")
        return None


def download_many(uniprot_ids, out_dir="alphafold_structures", version=4, format="pdb"):
    """
    Download AlphaFold structures for a list of UniProt IDs.
    """
    results = {"success": [], "not_found": [], "error": []}

    for uid in uniprot_ids:
        path = download_alphafold_structure(uid, out_dir=out_dir, version=version, format=format)
        if path:
            results["success"].append(uid)
        else:
            results["not_found"].append(uid)

    print(f"\nDone: {len(results['success'])} downloaded, "
          f"{len(results['not_found'])} not found")
    return results

# def get_alphafold_structure(uniprot_id, format="pdb"):
#     if format != "pdb" and format != "cif":
#         raise ValueError("format must be 'pdb' or 'cif'")
#     file_name_v4 = f"AF-{uniprot_id}-F1-model_v4.{format}.gz"
#     file_name_v6 = f"AF-{uniprot_id}-F1-model_v6.{format}.gz"
#     file_path_v4 = STRUCTURE_DIR / "alphafold" / file_name_v4
#     file_path_v6 = STRUCTURE_DIR / "alphafold" / file_name_v6
#     if file_path_v4.exists():
#         return file_path_v4
#     elif file_path_v6.exists():
#         return file_path_v6
#     else:
#         return None

# def uniprot_ids_with_saved_alphafold_structures(format="pdb"):
#     if format != "pdb" and format != "cif":
#         raise ValueError("format must be 'pdb' or 'cif'")
#     folder = STRUCTURE_DIR / "alphafold"
#     files = folder.glob(f"AF-*F1-model_v*.{format}.gz")
#     uniprot_ids = [file.stem.split("-")[1] for file in files]
#     uniprot_ids = list(set(uniprot_ids))
#     return uniprot_ids


if __name__ == "__main__":
    # single protein
    # download_alphafold_structure(ALL_UNIPROT_IDS[0], out_dir="structures")

    # # multiple proteins
    # download_many(ALL_UNIPROT_IDS, out_dir="structures")
    dir = STRUCTURE_DIR / "alphafold"
    not_found = []
    for uniprot_id in ALL_UNIPROT_IDS:
        file_path_v4 = dir / f"AF-{uniprot_id}-F1-model_v4.pdb.gz"
        file_path_v6 = dir / f"AF-{uniprot_id}-F1-model_v6.pdb.gz"
        if file_path_v4.exists() or file_path_v6.exists():
            continue
        else:
            not_found.append(uniprot_id)
            #download_alphafold_structure(uniprot_id, out_dir=dir, version=4, format="pdb")
            download_alphafold_structure(uniprot_id, out_dir=dir, version=6, format="pdb")
    print(f"{len(not_found)} AF structures not found")
    #print(not_found)