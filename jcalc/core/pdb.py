# Third Party imports
import math
import logging
from pathlib import Path
from Bio.PDB import (
    PDBParser,
    Selection,
    NeighborSearch,
    calc_dihedral
)


# Base imports
from jcalc.settings import HUGGINS_ELECTRO


class JCalcPdb:
    """ Class to store vicinal coupling constant (3JH,H) for a given PDB
        structure. Receives as input PDB filename and J input file with
        chosen JH,H
    """

    def __init__(self, pdb, j_input):

        parser = PDBParser()
        self.wkdir = Path.cwd()
        self.pdb = pdb.replace(".pdb", "")
        self.struct = parser.get_structure(pdb, pdb)
        self.atom_list = Selection.unfold_entities(self.struct, 'A')
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
        with open(str(self.j_input), "r") as file:
            for line in file:
                line = line.split("\t")
                line[-1] = line[-1].replace("\n", "")
                j_list.append(list(map(int, line)))
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
                atom_dict[atom.serial_number] = \
                 [atom.get_vector(), atom.element, atom.get_coord()]

        self.atom_dict = atom_dict

    def search_subs(
        self,
        hx_atom,
        cx_atom,
        cy_atom,
        hy_atom,
        j_atoms,
        chosen_j
    ):
        """ Description:

            Usage:

            Parameters:
        """

        subs_list = []
        # j_atoms = [j[0] for j in j_atoms]
        center = self.atom_dict[cx_atom][2]
        ns = NeighborSearch(self.atom_list)
        neighbors = ns.search(center, 1.7)
        for neigh_atom in neighbors:
            if neigh_atom.serial_number not in j_atoms:
                subs_list.append(neigh_atom)

        for subs in subs_list:
            self.j_dict[chosen_j]["substituents"][subs.serial_number] = {
                "SY": self.atom_dict[subs.serial_number][0],
                "HX": self.atom_dict[hx_atom][0],
                "CX": self.atom_dict[cx_atom][0],
                "CY": self.atom_dict[cy_atom][0],
                "element": self.atom_dict[subs.serial_number][1]
            }

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
            self.j_dict = j_dict

            # Search subs for CX and create subs dict
            self.search_subs(hx_atom=j[0],
                             cx_atom=j[1],
                             cy_atom=j[2],
                             hy_atom=j[3],
                             j_atoms=j,
                             chosen_j=chosen_j
                             )
            # Search subs for CY and create subs dict
            self.search_subs(hx_atom=j[3],
                             cx_atom=j[2],
                             cy_atom=j[1],
                             hy_atom=j[0],
                             j_atoms=j,
                             chosen_j=chosen_j
                             )

    def calc_subs_coupling(self, HX, CX, CY, SY, dihedral, element):
        """ Description:
              Given an JCalcPdb struct and atom vectors, calculate
              J pertubation from substitute atoms

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
        subs_dih = calc_dihedral(HX, CX, CY, SY)

        if subs_dih >= 0:
            return huggins_constant * \
                (0.56 + (-2.32 * (math.cos(math.radians((dihedral * -1) +
                 (17.9 * huggins_constant))) ** 2)))

        else:
            return huggins_constant * \
                (0.56 + (-2.32 * (math.cos(math.radians(dihedral +
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
            subs_value += \
                self.calc_subs_coupling(
                    HX=subs["HX"],
                    CX=subs["CX"],
                    CY=subs["CY"],
                    SY=subs["SY"],
                    dihedral=self.j_dict[chosen_j]["dih"],
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

    def write_pdb_results(self):
        """ Description:

            Usage:

            Parameters:
        """

        out_file = self.wkdir.joinpath(f"{self.pdb}_J_values.tsv")
        with open(str(out_file), "w") as j_file:
            for j, j_value in self.j_values.items():
                j_value = self.j_values[j]
                j_file.write(f"{j}\t{round(j_value,2)}\n")
        logging.info(f"Output file path: {str(out_file.resolve())}")
