
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
