import os
import numpy as np
import sys
import math
import subprocess
import statistics
from Bio.PDB import *
from jcalc.jcalcpdb import JCalcPdb
from tools.settings import HUGGINS_ELECTRO,GROMACS_VERSION

class JCalcMd:
    """
    """

    def __init__(self, xtc, tpr, residue, suffix, skip, j_input):
        self.xtc = xtc
        self.tpr = tpr
        self.residue = residue
        self.suffix = suffix
        self.skip = skip
        self.j_input = j_input

    def create_frames(self):
        """ Description:
              Given a JCalcMd struct, separate its XTC file in PDB files,
              skipping n frames chosen (self.skip)

            Usage:
              JCalcMd.create_frames()
        """
        os.mkdir(f"frames{self.suffix}/")
        subprocess.call(f"echo 2 2 | {GROMACS_VERSION} trjconv -s {self.tpr} \
                          -f {self.xtc} -sep -skip {self.skip} \
                          -o frames{self.suffix}/frame_.pdb -pbc mol -center",
                          shell=True
                       )
        frames = os.listdir(f"frames{self.suffix}/")
        frames = sorted(frames,
                        key=lambda x: int(x.split("_")[1].split(".")[0])
                       )
        self.frames = frames

    def add_hydrogen(self):
        """ Description:
              Given a JCalcMd struct, add hydrogens to all frames from
              the simulation

            Usage:
              JCalcMd.add_hydrogen()
        """

        for pdb in self.frames:
            subprocess.call(f"echo 3 | {GROMACS_VERSION} pdb2gmx -quiet \
                             -f frames{self.suffix}/{pdb} \
                             -o frames{self.suffix}/{pdb}_hydro.pdb \
                             -ff add_hydrogen", shell=True
                           )
            subprocess.call(f"rm *.top", shell=True)
            subprocess.call(f"rm *.itp", shell=True)

        frames = os.listdir(f"frames{self.suffix}/")
        hydro_frames = []
        for new in frames:
            if "hydro" in new:
                hydro_frames.append(new)

        frames = sorted(hydro_frames,
                        key=lambda x: int(x.split("_")[1].split(".")[0])
                       )
        self.frames = frames

    def calc_md_j(self):
        """ Description:

            Usage:
        """

        all_j_values = {}
        n_frames = 0

        for pdb in self.frames:
            j_struct = JCalcPdb(pdb=f"frames{self.suffix}/{pdb}",
                                j_input=self.j_input
                               )
            j_struct.get_atoms_vector()
            j_struct.create_j_dict()
            j_struct.calc_all_j()
            all_j_values[pdb] = j_struct
            n_frames += 1

        self.n_frames = n_frames
        self.all_j_values = all_j_values

        # Get J names
        j_names = []
        first_frame = list(self.all_j_values.keys())[0]
        for j in all_j_values[first_frame].j_values:
            j_names.append(j)
        self.j_names = j_names

    def calc_statistics(self):
        """ Description:
              Given a JCalcMd struct, calculate all statistics from J values
              calculated through Molecular Dynamics

            Usage:
              JCalcMd.calc_statistics()
        """

        statistics_dict = {}

        for j in self.j_names:
            statistics_dict[j] = []

        for pdb in self.all_j_values:
            for j_name,j_value in self.all_j_values[pdb].j_values.items():
                statistics_dict[j_name].append(j_value)

        # Now, calc statistics
        # Mean calc
        mean_results = {}
        for j in statistics_dict:
            mean_results[j] = statistics.fmean(statistics_dict[j])
        self.mean_results = mean_results

        # Stdev calc
        stdev_results = {}
        for j in statistics_dict:
            mean_results[j] = statistics.stdev(statistics_dict[j])

        self.stdev_results = stdev_results

    def write_statistics(self, out_name):
        """ Description:
              Given an JCalcMd and a output filename, returns statistical
              results from coupling constant (J) values

            Usage:
              JCalcMd.write_statistics("statistical_results.txt")

            Parameters:
              out_name:
                string, statitics output file name
        """

        with open(out_name, "w") as out:
            for j,mean_value in self.mean_results.items():
                out.write(f"{j}_mean:\t{mean_value}\n")
                out.write(f"{j}_stdev:\t{self.stdev_results[j]}\n")

    def write_j_values(self):
        """ Description:
              Given a JCalcMd, writes all J values for all frames in the
              Molecular Dynamics

            Usage:
            JCalcMd.write_j_values()
        """

        for j in self.j_names:
            with open(f"{j}_values.tsv","w") as j_file:
                for pdb in self.all_j_values:
                    j_value = self.all_j_values[pdb].j_values[j]
                    j_file.write(f"{pdb}\t{round(j_value,2)}\n")

    def write_results(self, stats_filename):

        """ Description:
              Given a JCalcMd struct, writes every output possible, being:
              statistics output;
              J values through frames;

            Usage:
              JCalcMd.write_results(stats_filename="statistical_results.txt")

            Parameters:
              stats_filename:
                string, output statistics file name
        """

        # Write statistics results
        self.write_statistics(out_name=stats_filename)

        # Write J values results through Molecular Dynamics
        self.write_j_values()
