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
import argparse



if __name__ == "main":
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xtc", "--xtc_filename", dest="xtc",
                        help=""".xtc filename from molecular dynamics you \
want to analyze (Needs to be in current directory)"""
                       )

    parser.add_argument("-t", "--tpr", "--tpr_filename", dest="tpr",
                        help=""".tpr filename from molecular dynamics you \
want to analyze (Needs to be in current directory)"""
                       )

    parser.add_argument("-r", "--res", "--residue_name", dest="residue",
                        help="""Residue name you want to calculate vicinal \
coupling constant (3JH,H) during molecular dynamics"""
                       )

    parser.add_argument("--suf", "--suffix", dest="suffix",
                        help="""Suffix of analysis, used when running multiple \
3JH,H calculations in same directory""", required=False
                       )

    parser.add_argument("--skip", dest="skip", type=int
                        help="""Number of frames you want to skip when calcu\
lating 3JH,H values thorugh molecular dynamics""", required=False,
                        default=100
                   )

    parser.add_argument("--j", "--j_input", dest="j_input", type=str
                        help="""J input filename to know which coupling \
constants will be calculated. Read documentation for information on format"""
                       )

    args = parser.parse_args()


    teste = JCalcMd(xtc=args.xtc,
                    tpr=args.tpr,
                    residue=args.residue,
                    suffix=args.suffix,
                    skip=args.skip,
                    j_input=args.j_input
                   )
    teste.create_frames()
    teste.add_hydrogen()
    teste.calc_md_j()
    teste.calc_statistics()
    teste.write_results()
