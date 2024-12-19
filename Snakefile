
from snakemake.utils import min_version
min_version("6.0")

import sys
import pathlib

sys.path.append("workflow/pypsa-earth")

rule copy_custom_powerplants:
    input: 
        "data/custom_powerplants.csv"
    output:
        "workflow/pypsa-earth/data/custom_powerplants.csv"
    shell:
        "cp {input} {output}"

rule retrieve_data:
    params:
        gdrive_url = "https://drive.google.com/drive/folders/1LwSoZDtnyUx5ki9SmBvdlGW3QwWjs4rA?usp=drive_link",
        output_path = "analysis/gdrive_data/"
    script:
        "analysis/scripts/download_from_gdrive.py"
