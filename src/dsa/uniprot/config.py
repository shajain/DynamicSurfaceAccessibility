from dsa.config import PROJECT_ROOT
from dsa.ms_dsa.ms_dsa import MS_DSA_Builder
import random

sequence_file = PROJECT_ROOT / "src" / "dsa" / "data" / "sequences.json"
ALL_UNIPROT_IDS = MS_DSA_Builder.get_all_uniprot_ids()
example_uniprot_id = ALL_UNIPROT_IDS[random.randint(0, len(ALL_UNIPROT_IDS)-1)]