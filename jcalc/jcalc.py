import os
import sys
import math
import subprocess
import statistics
from Bio.PDB import *
from jcalc.settings import *


class JCalcPDB:
    """
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
                JCalc.parse_j_list()
        """

        j_list = []
        with open(self.j_input, "r") as file:
            for line in file:
                line = line.split("\t")
                line[-1].replace("\n","")
                j_list.append(line)

        self.j_list = j_list

    def get_atoms_vector(self):
        """ Description:

            Usage:

            Parameters:
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

            Usage:

            Parameters:
        """

        j_dict = {}

        for j in self.j_list:
            chosen_j = f"{j[0],j[3]}"
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

    def calc_subs_coupling(self, HX, CX, CY, SY, dihedral, element):
        """ Description:

            Usage:

            Parameters:
        """

        huggins_const = huggins_electro[element]
        subs_dih = calc_dihedral(HX,CX,CY,SY)

        if subst_dih >= 0:
            return huggins_const * \
            (0.56  + (-2.32 * \
            (math.cos(math.radians((dihedral * -1) + \
            (17.9 * huggins_const))) ** 2)))

        else:
            return huggins_const * \
            (0.56  + (-2.32 * \
            (math.cos(math.radians(dihedral + \
            (17.9 * huggins_const))) ** 2)))

    def calc_j_h_h(self, chosen_j):
        """ Description:

            Usage:

            Parameters:
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
        print(subs_value)
        dih_radians = math.radians(chosen_j.dihedral)
        j_value = (13.86 * (math.cos(dih_radians)) ** 2) + \
                  (-0.81 * math.cos(dih_radians)) + subs_value
        return j_value


teste = JCalcPDB(pdb="teste.pdb", j_input="teste.inp")
teste.get_atoms_vector()
teste.create_j_dict()
print(teste.calc_j_h_h("H1,H2"))
