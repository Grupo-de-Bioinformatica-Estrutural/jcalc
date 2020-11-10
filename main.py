import os
import numpy as np
import sys
import math
import subprocess
import statistics
from Bio.PDB import *
from tools.settings import HUGGINS_ELECTRO,GROMACS_VERSION
from jcalc.jcalcpdb import JCalcPdb
from jcalc.jcalcmd import JCalcMd

teste = JCalcMd(xtc="sta_200_201ns.xtc",
                tpr="sta_prod_1ms.tpr",
                residue="STA",
                suffix=1,
                skip=100,
                j_input="j_input.txt"
               )
teste.create_frames()
teste.add_hydrogen()
teste.calc_md_j()
teste.calc_statistics()
teste.write_results()
