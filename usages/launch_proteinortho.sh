#!/bin/bash -i
conda activate mybase
perl ./proteinortho -project=projectname --cpus=20 --ram=250000 ../projectname/*faa

