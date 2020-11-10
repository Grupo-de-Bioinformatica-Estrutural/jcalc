import os
import numpy as np
import sys
import math
import subprocess
import statistics
from Bio.PDB import *
from jcalc.settings import HUGGINS_ELECTRO,GROMACS_VERSION


class JCalcPdb:
    """ Class to store vicinal coupling constant (3JH,H) for a given PDB
        structure. Receives as input PDB filename and J input file with
        chosen JH,H
    """

    def __init__(self, pdb, j_input):

        parser = PDBParser()
        self.struct = parser.get_structure(pdb, pdb)
        self.j_input = j_input
        self.parse_j_list()

    def parse_j_list(self):
        """ Description:
              Given a JCalc struct, create new attribute with all chosen
              Vicinal Coupling Constant to be calculated, being the attribute
              a list where each item is:
              list[0] = first proton name  ("H1")
              list[1] = second proton name ("H2")

            Usage:
                JCalcPdb.parse_j_list()
        """

        j_list = []
        n_j = 0
        with open(self.j_input, "r") as file:
            for line in file:
                line = line.split("\t")
                line[-1] = line[-1].replace("\n","")
                j_list.append(line)
                n_j += 1

        self.n_j = n_j
        self.j_list = j_list

    def get_atoms_vector(self):
        """ Description:
              Given an JCalcPdb struct, get all atoms vectors and elements
              from PDB struct

            Usage:
              JCalcPdb.get_atoms_vector()
        """

        structure = self.struct
        atom_dict = {}
        # Pegar so o residuo que quero como struct
        for residue in structure.get_residues():
            for atom in residue:
                atom_dict[atom.get_id()] = [atom.get_vector(), atom.element]
        self.atom_dict = atom_dict

    def create_j_dict(self):
        """ Description:
              Given an JCalcPdb struct, create J dictioary with chosen 3JH,H
              information to calc all J values

            Usage:
              JcalcPdb.create_j_dict()
        """

        j_dict = {}

        for j in self.j_list:
            chosen_j = f"{j[0]},{j[3]}"
            j_dict[chosen_j] = {}
            j_dict[chosen_j]["HX"] = self.atom_dict[j[0]][0]
            j_dict[chosen_j]["CX"] = self.atom_dict[j[1]][0]
            j_dict[chosen_j]["CY"] = self.atom_dict[j[2]][0]
            j_dict[chosen_j]["HY"] = self.atom_dict[j[3]][0]
            j_dict[chosen_j]["dih"] = calc_dihedral(self.atom_dict[j[0]][0],
                                                    self.atom_dict[j[1]][0],
                                                    self.atom_dict[j[2]][0],
                                                    self.atom_dict[j[3]][0]
                                                   )
            j_dict[chosen_j]["dih"] = math.degrees(j_dict[chosen_j]["dih"])
            j_dict[chosen_j]["substituents"] = {}
            j_dict[chosen_j]["substituents"][j[4]] = \
            {"SY": self.atom_dict[j[4]][0],
             "HX": self.atom_dict[j[3]][0],
             "CX": self.atom_dict[j[2]][0],
             "CY": self.atom_dict[j[1]][0],
             "element": self.atom_dict[j[4]][1]
            }
            j_dict[chosen_j]["substituents"][j[5]] = \
            {"SY": self.atom_dict[j[5]][0],
             "HX": self.atom_dict[j[3]][0],
             "CX": self.atom_dict[j[2]][0],
             "CY": self.atom_dict[j[1]][0],
             "element": self.atom_dict[j[5]][1]
            }
            j_dict[chosen_j]["substituents"][j[6]] = \
            {"SY": self.atom_dict[j[6]][0],
             "HX": self.atom_dict[j[0]][0],
             "CX": self.atom_dict[j[1]][0],
             "CY": self.atom_dict[j[2]][0],
             "element": self.atom_dict[j[6]][1]
            }
            j_dict[chosen_j]["substituents"][j[7]] = \
            {"SY": self.atom_dict[j[7]][0],
             "HX": self.atom_dict[j[0]][0],
             "CX": self.atom_dict[j[1]][0],
             "CY": self.atom_dict[j[2]][0],
             "element": self.atom_dict[j[7]][1]
            }

        self.j_dict = j_dict

    def calc_subs_coupling(self, HX, CX, CY, SY, dihedral, element):
        """ Description:
              Given an JCalcPdb struct and atom vectors, calculate J pertubation
              from substitute atoms

            Usage:
              JcalcPdb.calc_subs_coupling(vector_HX, vector_CX, vector_CY,
                                          vector_SY, HX-CX-CY-XY_dih,
                                          substitue_element
                                         )

            Parameters:
              HX:
                Atom coordinates vector, HX from HX-CX-CY-SY fragment
              CX:
                Atom coordinates vector, CX from HX-CX-CY-SY fragment
              CY:
                Atom coordinates vector, CY from HX-CX-CY-SY fragment
              SY:
                Atom coordinates vector, SY from HX-CX-CY-SY fragment
              dihedral:
                float, dihedral from JH,H fragment HX-CX-CY-XY in radians
              element:
                string, atom element from substitute atom (H, O, N)
        """

        huggins_constant = HUGGINS_ELECTRO[element]
        subs_dih = calc_dihedral(HX,CX,CY,SY)

        if subs_dih >= 0:
            return huggins_constant * \
            (0.56  + (-2.32 * \
            (math.cos(math.radians((dihedral * -1) + \
            (17.9 * huggins_constant))) ** 2)))

        else:
            return huggins_constant * \
            (0.56  + (-2.32 * \
            (math.cos(math.radians(dihedral + \
            (17.9 * huggins_constant))) ** 2)))

    def calc_j_h_h(self, chosen_j):
        """ Description:
              Given an JCalcPdb struct and a chosen 3JH,H, calculate it's J
              value and return it

            Usage:
              JCalcPdb.calc_j_h_h(chosen_j="H1,H2")

            Parameters:
              chosen_j:
                string, chosen 3JH,H value from input couplings contants
                example: "H1,H2"
        """

        subs_value = 0
        for subs in self.j_dict[chosen_j]["substituents"].values():
            subs_value += self.calc_subs_coupling(HX=subs["HX"],
                                                  CX=subs["CX"],
                                                  CY=subs["CY"],
                                                  SY=subs["SY"],
                                                  dihedral=\
                                                  self.j_dict[chosen_j]["dih"],
                                                  element=subs["element"]
                                                 )

        dih_radians = math.radians(self.j_dict[chosen_j]["dih"])
        j_value = (13.86 * (math.cos(dih_radians)) ** 2) + \
                  (-0.81 * math.cos(dih_radians)) + subs_value
        return j_value

    def calc_all_j(self):
        """ Description:
              Given an JcalcPdb struct, call all 3JH,H values given as inputs
              from j_input attribute

            Usage:
              JcalcPdb.calc_all_j()
        """

        j_values = {}
        for j in self.j_dict.keys():
            j_values[j] = self.calc_j_h_h(j)

        self.j_values = j_values


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
                statistics_values[j_name].append(j_value)

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
                    j_file.write(f"{pdb}\t{j_value}\n")

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
