from pathlib import Path

def get_project_root() -> Path:
    """
    Returns absolute path to PATHFINDER project root
    :return:
    """
    pack_root = Path(__file__).resolve().parent
    return pack_root.parent.parent

PROJECT_ROOT = get_project_root()
SOURCE_ROOT = PROJECT_ROOT / "src"
PATH_ROOT = PROJECT_ROOT / "src" / "prism"

# =====================
# 1st level directories
# =====================
DATABASE = PATH_ROOT / "database"
GENE_DATABANK = PATH_ROOT / "gene_databank"
REGEN_CANDIDATES = PATH_ROOT / "ReGen_candidates"
SCREEN_RESULTS = PATH_ROOT / "Screen_results"
DATASIFT_CONFIGS = PATH_ROOT / "DataSift_configs"
MODEL_FOLDER = PATH_ROOT / "model_folder"

# =====================
# embedded directories
# =====================
PWM_DATABASE = DATABASE / "pwm_database"
DNA_MOTIFS = PWM_DATABASE / "DNA_motifs"
AA_MOTIFS = PWM_DATABASE / "AA_motifs"

# blosum matrix
BLOSUM62 = DATABASE / "blosum62.mat"

# config location
CONFIG = PATH_ROOT / "config.yaml"





