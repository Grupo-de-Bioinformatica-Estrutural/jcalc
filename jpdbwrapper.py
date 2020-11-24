#!/usr/bin/env python3
__author__ = "Joao Luiz de Meirelles"

import os
import numpy as np
import sys
import math
import subprocess
import statistics
from Bio.PDB import *
from tools.settings import HUGGINS_ELECTRO,GROMACS_VERSION
from jcalc.jcalcpdb import JCalcPdb
from jcalc.jcalcpdb import main
from jcalc.jcalcmd import JCalcMd
import argparse

if __name__ == "__main__":
    main()
