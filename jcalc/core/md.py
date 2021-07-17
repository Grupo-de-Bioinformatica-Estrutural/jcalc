# Third Party imports
import os
import subprocess
import statistics
import logging
from pathlib import Path
import multiprocessing


# Base imports
from jcalc.core.pdb import JCalcPdb
from jcalc.settings import GROMACS_VERSION


class JCalcMd:
    """
    """

    def __init__(self, xtc, tpr, suffix, skip, j_input):

        self.wkdir = Path.cwd()
        self.xtc = xtc
        self.tpr = tpr
        self.suffix = suffix
        self.skip = skip
        self.j_input = self.wkdir.joinpath(j_input)
        self.frames_dir = self.wkdir.joinpath(f"frames{self.suffix}")

    def create_frames(self):
        """ Description:
              Given a JCalcMd struct, separate its XTC file in PDB files,
              skipping n frames chosen (self.skip)

            Usage:
              JCalcMd.create_frames()
        """
        frames_dir = self.frames_dir
        if frames_dir.exists():
            logging.error(f"Dir {str(frames_dir.resolve())} exists, \
please rename it or remove it")
        else:
            frames_dir.mkdir()
        logging.info(f"Creating frames, path: {str(frames_dir.resolve())}")

        subprocess.call(f"echo non-Water non-Water | \
                          {GROMACS_VERSION} trjconv -s {self.tpr} \
                          -f {self.xtc} -sep -skip {self.skip} \
                          -o {str(frames_dir.resolve())}/frame_.pdb -pbc mol \
                          -center",
                        shell=True
                        )
        frames = os.listdir(str(frames_dir.resolve()))
        frames = sorted(frames,
                        key=lambda x: int(x.split("_")[1].split(".")[0])
                        )
        self.frames = frames

    def rename_hydro(self, pdb):
        """ Description:
              Given an PDB file, rename all added hydrogen by obabel so
              BioPDBParser can accept its format

            Usage:
              JCalcMd.rename_hydro()

            Parameters:
              pdb:
                string, pdb filename
        """

        old_h = []
        new_file_lines = []
        with open(pdb, "r") as file:
            for line in file:
                cur_line = line.split()
                try:
                    if "H" in cur_line[2] and \
                     (cur_line[0] == "ATOM" or cur_line[0] == "HETATM"):
                        if cur_line[2] != "H" and cur_line[2] not in old_h:
                            old_h.append(cur_line[2])

                    new_file_lines.append(line)
                except Exception as error:
                    new_file_lines.append(line)
                    logging.warning(f"Error {error} at renaming hydrogen")
                    continue

        # Now, we have all hidrogens names, and we can add new ones without
        # messing with original hydrogens
        counter = 1
        with open(pdb, "w") as file:
            for line in new_file_lines:
                cur_line = line.split()
                try:
                    if cur_line[2] == "H" and \
                     (cur_line[0] == "ATOM" or cur_line[0] == "HETATM"):
                        h_name = f"H{counter}"
                        while h_name in old_h:
                            counter += 1
                            h_name = f"H{counter}"
                        line = line[:12] + f"{h_name:^4s}" + f" {line[17:]}"
                        old_h.append(h_name)
                        file.write(line)
                    else:
                        file.write(line)
                except Exception as error:
                    logging.warning(f"Error {error} at renaming hydrogen")
                    file.write(line)

    def add_hydrogen(self):
        """ Description:
              Given a JCalcMd struct, add hydrogens to all frames from
              the simulation

            Usage:
              JCalcMd.add_hydrogen()
        """
        jobs = []
        logging.info("Adding hydrogen to frames (GROMOS non-polar hydrogen)")
        for pdb in self.frames:
            process = multiprocessing.Process(
                target=self.cmd_addhydro,
                args=(pdb,)
                )
            jobs.append(process)
            process.start()
        for j in jobs:
            j.join()

    def cmd_addhydro(self, pdb):
        """ Description:

            Usage:
        """
        frames_dir = self.frames_dir
        pdb_file = f"{str(frames_dir.resolve())}/{pdb}"
        cmd_add = f"obabel -ipdb {pdb_file} \
-opdb -O {pdb_file} -h"
        subprocess.call(cmd_add, shell=True)
        self.rename_hydro(f"{pdb_file}")

    def calc_md_j(self):
        """ Description:

            Usage:
        """

        logging.info("Calculating J values from MD")
        frames_dir = self.frames_dir
        all_j_values = {}
        n_frames = 0

        for pdb in self.frames:
            j_struct = JCalcPdb(pdb=f"{str(frames_dir.resolve())}/{pdb}",
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

        logging.info("Calculating Statistics from J values")
        statistics_dict = {}

        for j in self.j_names:
            statistics_dict[j] = []

        for pdb in self.all_j_values:
            for j_name, j_value in self.all_j_values[pdb].j_values.items():
                statistics_dict[j_name].append(j_value)

        # Now, calc statistics
        # Mean calc
        mean_results = {}
        for j in statistics_dict:
            mean_results[j] = statistics.mean(statistics_dict[j])
        self.mean_results = mean_results

        # Stdev calc
        stdev_results = {}
        for j in statistics_dict:
            stdev_results[j] = statistics.stdev(statistics_dict[j])

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

        with open(str(out_name), "w") as out:
            for j, mean_value in self.mean_results.items():
                out.write(f"{j}_mean:\t{round(mean_value,2)}\n")
                out.write(f"{j}_stdev:\t{round(self.stdev_results[j],2)}\n")

    def write_j_values(self):
        """ Description:
              Given a JCalcMd, writes all J values for all frames in the
              Molecular Dynamics

            Usage:
            JCalcMd.write_j_values()
        """

        for j in self.j_names:
            out_file = self.wkdir.joinpath(f"{j}_values.tsv")
            logging.info(f"Writing J {j} values through MD, path: \
{str(out_file)}")
            with open(str(out_file), "w") as j_file:
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
        stats_file = self.wkdir.joinpath(f"{stats_filename}")
        self.write_statistics(out_name=stats_file)
        logging.info(f"Writing J statistics resuls, path: \
{str(stats_file.resolve())}")

        # Write J values results through Molecular Dynamics
        self.write_j_values()
