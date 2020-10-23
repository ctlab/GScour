#!/bin/bash -i
conda activate mybase
perl ./proteinortho -project=acfhp_synteny --cpus=32 --ram=250000 --synteny ../acfhp_synteny/gfffiles/*faa

