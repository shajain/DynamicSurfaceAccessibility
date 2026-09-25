import gzip
import requests
import xml.etree.ElementTree as ET
import pandas as pd
from pathlib import Path
from dsa.misc.downloader import Downloader
from dsa.misc.store import CompressedPickleStore, Store
from dsa.binding_interfaces.config import SIFT_alignment_file, PDBE_ALIGNMENT_DIR, ALIGNMENT_DIR
from dsa.binding_interfaces.mappings.uniprot_mappings import UniprotStructureMappingsStore
PDB_STRUCTURES_TO_PROCESS = UniprotStructureMappingsStore("pdb").get_mapper().get_mapped_structures()
from dsa.binding_interfaces.factory import ClassFactory
from dsa.binding_interfaces.structure.structure_store import StructureFileStore
from dsa.misc.compute import ComputeManager
import pandas as pd
from dsa.misc.store import SingleFileTabularStore

class SIFTMappingsStore(SingleFileTabularStore):
    def __init__(self):
        super().__init__(directory=ALIGNMENT_DIR, file_name=SIFT_alignment_file, dtype='str')
        #self.downloader = PDBeDownloader()

    def load_mappings(self, structure_id: str) -> pd.DataFrame:
        return self.data[self.data["PDB"] == structure_id]
        


class PDBeUniprotAlignmentStore(CompressedPickleStore):
    def __init__(self):
        super().__init__(directory=PDBE_ALIGNMENT_DIR, extension="alignment", object_type=pd.DataFrame)
        #self.downloader = PDBeDownloader()

    def load_mappings(self, structure_id: str) -> pd.DataFrame:
        return self.load(structure_id)

    # def download(self, structure_id: str) -> pd.DataFrame:
    #     data = self.downloader.download(structure_id)
    #     if data is None:
    #         return None
    #     self.dump(structure_id, data)
    #     return data


class PDBeDownloader(Downloader):
    def __init__(self):
        base_url = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{structure_id}.xml.gz"
        super().__init__(base_url=base_url)

    def build_url(self, structure_id: str) -> str:
        return self.base_url.format(structure_id=structure_id.lower())

    def process_response(self, response, structure_id: str) -> pd.DataFrame:
        try:
            xml_content = gzip.decompress(response.content)
        except Exception as e:
            print(f"Error decompressing {structure_id}: {e}")
            return None
        root = ET.fromstring(xml_content)
        ns = {'sifts': 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'}
        residue_mappings = []

        for entity in root.findall('.//sifts:entity', ns):
            entity_id = entity.get('entityId')  # NOTE: alphabetical index only, NOT the real chain ID — kept for debugging/traceability only

            for residue in entity.findall('.//sifts:residue', ns):
                pdb_res = residue.find('sifts:crossRefDb[@dbSource="PDB"]', ns)
                # Use findall, not find — some residues can carry more than one UniProt xref (e.g. canonical + isoform context)
                unp_res_list = residue.findall('sifts:crossRefDb[@dbSource="UniProt"]', ns)

                if pdb_res is None or not unp_res_list:
                    continue

                pdb_chain_id = pdb_res.get('dbChainId')  # this is the REAL author chain ID (e.g. "A", "B")

                for unp_res in unp_res_list:
                    uniprot_accession = unp_res.get('dbAccessionId')  # e.g. "P54748" or "P54748-3" for isoform 3
                    is_isoform = bool(uniprot_accession) and '-' in uniprot_accession
                    canonical_accession = uniprot_accession.split('-')[0] if uniprot_accession else None
                    isoform_number = (
                        uniprot_accession.split('-')[1] if is_isoform else None
                    )

                    residue_mappings.append({
                        'pdbe_record_id':        entity_id,          # kept for reference only, don't rely on this as chain ID
                        'pdb_chain_id':          pdb_chain_id,       # real author chain ID
                        'pdb_resnum':            pdb_res.get('dbResNum'),
                        'pdb_resname':           pdb_res.get('dbResName'),
                        'uniprot_id':            canonical_accession, # canonical accession
                        'uniprot_isoform_id':    uniprot_accession,  # full accession, isoform-suffixed if applicable,
                        'is_isoform':            is_isoform,
                        'isoform_number':        isoform_number,
                        'unp_resnum':            unp_res.get('dbResNum'),
                        'unp_resname':           unp_res.get('dbResName'),
                    })

        if residue_mappings:
            residue_mappings = pd.DataFrame(residue_mappings)
        else:
            residue_mappings = None
        return residue_mappings



def main(update_interval = 25, show_progress = True):
    database = "pdb"
    structure_ids = PDB_STRUCTURES_TO_PROCESS
    structure_file_store = ClassFactory.create(StructureFileStore, database)
    structure_ids = structure_file_store.get_saved_structures(from_structures=structure_ids)
    
    pdbe_uniprot_alignment_store = PDBeUniprotAlignmentStore()

    download_pdbe_alignment = lambda structure_id: PDBeDownloader().download(structure_id)

    compute_manager = ComputeManager(structure_ids, 
                                     compute_function=download_pdbe_alignment, 
                                     store=pdbe_uniprot_alignment_store,
                                     update_interval=update_interval,
                                     show_progress=show_progress)
    compute_manager.compute()

        

if __name__ == "__main__":
    main(update_interval=25, show_progress=True)



# PDB_UNIPROT_MAP_DIR = Path(__file__).parent / "pdbe_uniprot_mappings"
# PDB_UNIPROT_MAP_DIR.mkdir(parents=True, exist_ok=True)

# def download_mappings(structure_id: str):
#     url = f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{structure_id.lower()}.xml.gz"
#     response = requests.get(url)
#     try:
#         xml_content = gzip.decompress(response.content)
#     except Exception as e:
#         print(f"Error decompressing {structure_id}: {e}")
#         return None

#     root = ET.fromstring(xml_content)
#     ns = {'sifts': 'http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd'}
#     residue_mappings = []

#     for entity in root.findall('.//sifts:entity', ns):
#         entity_id = entity.get('entityId')  # NOTE: alphabetical index only, NOT the real chain ID — kept for debugging/traceability only

#         for residue in entity.findall('.//sifts:residue', ns):
#             pdb_res = residue.find('sifts:crossRefDb[@dbSource="PDB"]', ns)
#             # Use findall, not find — some residues can carry more than one UniProt xref (e.g. canonical + isoform context)
#             unp_res_list = residue.findall('sifts:crossRefDb[@dbSource="UniProt"]', ns)

#             if pdb_res is None or not unp_res_list:
#                 continue

#             pdb_chain_id = pdb_res.get('dbChainId')  # this is the REAL author chain ID (e.g. "A", "B")

#             for unp_res in unp_res_list:
#                 uniprot_accession = unp_res.get('dbAccessionId')  # e.g. "P54748" or "P54748-3" for isoform 3
#                 is_isoform = bool(uniprot_accession) and '-' in uniprot_accession
#                 canonical_accession = uniprot_accession.split('-')[0] if uniprot_accession else None
#                 isoform_number = (
#                     uniprot_accession.split('-')[1] if is_isoform else None
#                 )

#                 residue_mappings.append({
#                     'pdbe_entity_id':        entity_id,          # kept for reference only, don't rely on this as chain ID
#                     'pdb_chain_id':          pdb_chain_id,       # real author chain ID
#                     'pdb_resnum':            pdb_res.get('dbResNum'),
#                     'pdb_resname':           pdb_res.get('dbResName'),
#                     'uniprot_id':            canonical_accession, # canonical accession
#                     'uniprot_isoform_id':    uniprot_accession,  # full accession, isoform-suffixed if applicable,
#                     'is_isoform':            is_isoform,
#                     'isoform_number':        isoform_number,
#                     'unp_resnum':            unp_res.get('dbResNum'),
#                     'unp_resname':           unp_res.get('dbResName'),
#                 })

#     if residue_mappings:
#         pd.DataFrame(residue_mappings).to_csv(
#             PDB_UNIPROT_MAP_DIR / f"{structure_id}.tsv", sep="\t", index=False
#         )
#     return residue_mappings