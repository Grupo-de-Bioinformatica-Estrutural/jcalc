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
from jcalc.jcalcmd import JCalcMd
import argparse
import logging
# Log configurations
logging.basicConfig(filename="jcalc.log", filemode="w",
                    format="%(asctime)s - %(message)s",
                    datefmt='%m/%d/%Y %H:%M:%S',
                    level=logging.INFO
                   )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--x", "--xtc", "--xtc_filename", dest="xtc",
                        help=""".xtc filename from molecular dynamics you \
want to analyze (Needs to be in current directory)""", required=False
                       )

    parser.add_argument("--t", "--tpr", "--tpr_filename", dest="tpr",
                        help=""".tpr filename from molecular dynamics you \
want to analyze (Needs to be in current directory)""", required=False
                       )

    parser.add_argument("--r", "--res", "--residue_name", dest="residue",
                        help="""Residue name you want to calculate vicinal \
coupling constant (3JH,H) during molecular dynamics""", required=False
                       )

    parser.add_argument("--suf", "--suffix", dest="suffix",
                        help="""Suffix of analysis, used when running multiple \
3JH,H calculations in same directory""", required=False
                       )

    parser.add_argument("--skip", dest="skip", type=int,
                        help="""Number of frames you want to skip when calcu\
lating 3JH,H values thorugh molecular dynamics""", required=False,
                        default=100
                   )

    parser.add_argument("--j", "--j_input", dest="j_input", type=str,
                        help="""J input filename to know which coupling \
constants will be calculated. Read documentation for information on format""",
                        required=False
                       )

    parser.add_argument("--ff", "--ff_hydro", dest="ff_hydro", type=str,
                        help=""".ff dir with hydrogen information when adding \
hydrogen to implicit GROMOS simulations""", required=False
                       )


    parser.add_argument("--p", "--pdb", "--pdb_filename", dest="pdb",
                        help=""".pdb filename you want to calculate J values \
(Needs to be in current directory)""", required=False
                       )

    args = parser.parse_args()

    # If args.xtc exists, MD analysis is chosen
    if args.xtc:
        logging.info("Starting JCalc")
        logging.info(f"XTC chosen file: {args.xtc}")
        logging.info(f"TPR chosen file: {args.tpr}")
        logging.info(f"Chosen residue: {args.residue}")
        logging.info(f"Skip chosen: {args.skip}")
        logging.info(f"J input file: {args.j_input}")
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
        teste.write_results("statistical_results.txt")

    # Else, static structure analysis is chosen
    else:
        pdb_J = JCalcPdb(pdb=args.pdb, j_input=args.j_input)
        pdb_J.get_atoms_vector()
        pdb_J.create_j_dict()
        pdb_J.calc_all_j()
        pdb_J.write_pdb_results("teste.txt")
