import logging
import os.path
import re
from functools import cached_property
from glob import glob

import numpy as np
from ase import units
from matplotlib import pyplot as plt

from chemsmart.io.gaussian.output import Gaussian16Output
from chemsmart.io.orca.output import ORCAOutput
from chemsmart.utils.mixins import BaseFolder
from chemsmart.utils.repattern import (
    gaussian_dias_filename_point_with_fragment1,
    gaussian_dias_filename_point_with_fragment2,
    gaussian_dias_filename_point_without_fragment,
    gaussian_dias_filename_with_reactant,
    orca_dias_filename_point_with_fragment1,
    orca_dias_filename_point_with_fragment2,
    orca_dias_filename_point_without_fragment,
    orca_dias_filename_with_reactant,
)

logger = logging.getLogger(__name__)


class DIASOutputFolder(BaseFolder):
    """Analyse DI-AS model output folder.

    atom1: atom number 1 for plotting bond distance
    atom2: atom number 2 for plotting bond distance
    zero: bool to turn on/off reference to the lowest total energy as zero
    ref_file: user supplied reference file to act as minimum (energy in a.u. as from calculation).
    """

    def __init__(
        self,
        folder,
        atom1,
        atom2,
        zero=False,
        ref_file=None,
        outputname="dias",
    ):
        super().__init__(folder=folder)
        self.atom1 = atom1
        self.atom2 = atom2
        self.zero = zero  # zero the data with respect to relative total energy
        if zero:
            self.outputname = f"{outputname}_zero_ref"
        else:
            self.outputname = f"{outputname}_nonzero_ref"
        if ref_file is not None:
            self.outputname += "_user"
        self.ref_file = ref_file

        self.list_rc = self.get_reaction_coordinates()
        self.list_rel_total_energies = self.get_relative_total_energies()
        self.list_total_strain_energies = self.get_strain_energies()
        self.list_total_interaction_energies = self.get_interaction_energies()

    @property
    def irc_points(self):
        assert (
            len(self.full_molecule_energies_list)
            == len(self.fragment1_energies_list)
            == len(self.fragment2_energies_list)
        ), (
            f"Number of points along the IRC for full molecule {len(self.full_molecule_energies_list)} "
            f"should be the same as number of points for their corresponding fragments, f1, "
            f"{len(self.fragment1_energies_list)} and f2, {len(self.fragment2_energies_list)}"
        )
        return len(self.full_molecule_energies_list)

    @cached_property
    def reactant_energies_list(self):
        reactant_files_list = self._get_all_files_reactants()
        return self._get_all_energies(list_of_files=reactant_files_list)

    @cached_property
    def full_molecule_energies_list(self):
        full_molecule_files_list = self._get_all_files_full_molecule()
        return self._get_all_energies(list_of_files=full_molecule_files_list)

    @cached_property
    def fragment1_energies_list(self):
        fragment1_files_list = self._get_all_files_fragment1()
        logger.info(f"Fragment 1 files: {fragment1_files_list}")
        return self._get_all_energies(list_of_files=fragment1_files_list)

    @cached_property
    def fragment2_energies_list(self):
        fragment2_files_list = self._get_all_files_fragment2()
        return self._get_all_energies(list_of_files=fragment2_files_list)

    @property
    def ref_file_rel_energy(self):
        """Relative total energy of the reference file."""
        # convert reference file to a list
        if self.ref_file is not None:
            ref_file = [self.ref_file]
            ref_file_energy = self._get_all_energies(ref_file)
            assert (
                len(ref_file_energy) == 1
            ), "Only one reference file is expected!"
            return ref_file_energy[0] - sum(self.reactant_energies_list)
        return None

    @property
    def ref_file_rc(self):
        """Reaction coordinate of the reference file."""
        if self.ref_file is not None:
            ref_file = [self.ref_file]
            ref_file_rc = self._get_reaction_coordinates(ref_file)
            assert (
                len(ref_file_rc) == 1
            ), "Only one reference file is expected!"
            return ref_file_rc[0]
        return None

    def _get_all_files_full_molecule(self):
        """Subclass can implement."""
        raise NotImplementedError

    def _get_all_files_fragment1(self):
        """Subclass can implement."""
        raise NotImplementedError

    def _get_all_files_fragment2(self):
        """Subclass can implement."""
        raise NotImplementedError

    def _get_all_files_reactants(self):
        """Subclass can implement."""
        raise NotImplementedError

    def _get_all_energies(self, list_of_files):
        """Get all energies in kcal/mol as a list."""
        raise NotImplementedError

    def _get_ref_file(self):
        pass

    def get_relative_total_energies(self):
        """Total energy relative to the sum of reactants."""
        return [
            self.full_molecule_energies_list[i]
            - sum(self.reactant_energies_list)
            for i in range(self.irc_points)
        ]

    def get_strain_energies(self):
        # r1 should coincide with f1; r2 with f2
        return [
            (self.fragment1_energies_list[i] - self.reactant_energies_list[0])
            + (
                self.fragment2_energies_list[i]
                - self.reactant_energies_list[1]
            )
            for i in range(self.irc_points)
        ]

    def get_strain_energies_for_fragments(self):
        return [
            self.fragment1_energies_list[i] - self.reactant_energies_list[0]
            for i in range(self.irc_points)
        ], [
            self.fragment2_energies_list[i] - self.reactant_energies_list[1]
            for i in range(self.irc_points)
        ]

    def get_interaction_energies(self):
        # r1 should coincide with f1; r2 with f2
        return [
            self.full_molecule_energies_list[i]
            - self.fragment1_energies_list[i]
            - self.fragment2_energies_list[i]
            for i in range(self.irc_points)
        ]

    def get_reaction_coordinates(self):
        all_files_full_molecule_sorted = self._get_all_files_full_molecule()
        return self._get_reaction_coordinates(all_files_full_molecule_sorted)

    def _get_reaction_coordinates(self, list_of_filenames):
        """Computes distance between two atoms in each frame along IRC point, for a given list of filenames."""
        list_rc = []
        atom1 = (
            self.atom1 - 1
        )  # convert atom number from 1-index in gview to 0-index in python
        atom2 = self.atom2 - 1

        # all_files_full_molecule_sorted = self._get_all_files_full_molecule()
        for file in list_of_filenames:
            if file.endswith(".log"):
                gout = Gaussian16Output(filename=file)
                molecule = gout.molecule
            elif file.endswith(".out"):
                oout = ORCAOutput(filename=file)
                molecule = oout.molecule
            else:
                raise ValueError(
                    f"File {file} has unknown format. Acceptable formats are Gaussian .log or ORCA .out."
                )

            coordinates = molecule.positions
            symbols = molecule.chemical_symbols
            logger.info(
                f"Getting distance between atoms {symbols[atom1]}{self.atom1} and {symbols[atom2]}{self.atom2}"
            )
            distance = np.linalg.norm(coordinates[atom1] - coordinates[atom2])
            list_rc.append(distance)
        return list_rc

    def get_data(self):
        if self.zero:  # zeroing
            if self.ref_file is None:  # noqa: SIM108
                # no ref file given, take min of rel tot energy as zero
                lowest_rel_tot_e = min(self.list_rel_total_energies)
            else:
                # ref file is given, zero ref wrt to the ref file given
                lowest_rel_tot_e = self.ref_file_rel_energy
        else:
            # no zeroing reference, plot as given (whether ref file is given or not)
            lowest_rel_tot_e = 0.0
        rel_total_energies = [
            i - lowest_rel_tot_e for i in self.list_rel_total_energies
        ]
        total_strain_energies = [
            i - lowest_rel_tot_e for i in self.list_total_strain_energies
        ]
        total_interaction_energies = [
            i - lowest_rel_tot_e for i in self.list_total_interaction_energies
        ]
        return (
            rel_total_energies,
            total_strain_energies,
            total_interaction_energies,
        )

    def write_data(self):
        outputfile = f"{self.job_basename}_{self.outputname}_data.txt"
        outfile_path = os.path.join(self.folder, outputfile)

        (
            rel_total_energies,
            total_strain_energies,
            total_interaction_energies,
        ) = self.get_data()

        with open(outfile_path, "w") as f:
            f.write(
                "#   Reaction_coordinate    total    distortion    interaction\n"
            )
            for i in range(len(self.list_rc)):
                f.write(
                    f"{self.list_rc[i]:10.5f} {rel_total_energies[i]:10.5f} {total_strain_energies[i]:10.5f} "
                    f"{total_interaction_energies[i]:10.5f}\n"
                )

    def plot_dias(self, extrapolate=True, reversed=True, new_length=1000, k=3):
        from chemsmart.utils.utils import spline_data

        list_rc = self.list_rc
        (
            rel_total_energies,
            total_strain_energies,
            total_interaction_energies,
        ) = self.get_data()

        assert (
            len(list_rc)
            == len(rel_total_energies)
            == len(total_strain_energies)
            == len(total_interaction_energies)
        ), "Number of all points for full molecule and fragments must be the same!"
        plt.figure()
        font = {"family": "Arial", "weight": "bold"}
        plt.rc("font", **font)
        plt.rc("text")
        # matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

        # raw data points
        plt.plot(list_rc, rel_total_energies, "ko")
        plt.plot(list_rc, total_strain_energies, "bo")
        plt.plot(
            list_rc,
            total_interaction_energies,
            color="#7922BA",
            marker="o",
            linestyle="",
        )

        # also plot the reference file point
        if self.ref_file is not None:
            plt.plot(self.ref_file_rc, self.ref_file_rel_energy, "ko")

        # interpolate data
        if extrapolate:
            x, y_tot = spline_data(
                list_rc, rel_total_energies, new_length=new_length, k=k
            )
            x, y_strain = spline_data(
                list_rc, total_strain_energies, new_length=new_length, k=k
            )
            x, y_interaction = spline_data(
                list_rc, total_interaction_energies, new_length=new_length, k=k
            )

            plt.plot(x, y_tot, "k-", label=r"total")
            plt.plot(x, y_strain, "b-", label=r"distortion")
            plt.plot(
                x,
                y_interaction,
                color="#7922BA",
                linestyle="-",
                label=r"interaction",
            )

            min_y = min(y_interaction)
            max_y = max(y_strain)

        # no interpolation
        else:
            plt.plot(list_rc, rel_total_energies, "k-", label=r"total")
            plt.plot(list_rc, total_strain_energies, "b-", label=r"distortion")
            plt.plot(
                list_rc,
                total_interaction_energies,
                color="#7922BA",
                linestyle="-",
                label=r"interaction",
            )

            min_y = min(total_interaction_energies)
            max_y = max(total_strain_energies)

        x_range = max(list_rc) - min(list_rc)
        if reversed:
            plt.xlim(
                max(list_rc) + 0.1 * x_range, min(list_rc) - 0.1 * x_range
            )
        else:
            plt.xlim(
                min(list_rc) - 0.1 * x_range, max(list_rc) + 0.1 * x_range
            )

        y_range = max_y - min_y
        plt.ylim(min_y - 0.1 * y_range, max_y + 0.1 * y_range)

        plt.ylabel(
            r"$E$ / kcal mol$^{-1}$",
            color="black",
            fontsize=12,
            fontweight="bold",
        )
        plt.xlabel(
            r"Reaction coordinate / $\AA$",
            color="black",
            fontsize=12,
            fontweight="bold",
        )

        ax = plt.gca()

        # frame visibility
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(True)
        ax.spines["left"].set_visible(True)

        ax.yaxis.set_tick_params(length=5, direction="in")
        # legend box centre top
        ax.legend(
            loc="upper center",
            ncol=3,
            frameon=True,
            columnspacing=1.5,
            prop={"size": 12},
            handletextpad=0.3,
        )

        outputfile = f"{self.job_basename}_{self.outputname}_plot.pdf"
        outfile_path = os.path.join(self.folder, outputfile)
        plt.savefig(outfile_path, format="pdf", dpi=500)
        plt.show()


class GaussianDIASLogFolder(DIASOutputFolder):
    def __init__(
        self,
        folder,
        atom1,
        atom2,
        zero=False,
        ref_file=None,
        outputname="dias",
    ):
        super().__init__(
            folder=folder,
            atom1=atom1,
            atom2=atom2,
            zero=zero,
            ref_file=ref_file,
            outputname=outputname,
        )

    @property
    def job_basename(self):
        full_molecule_a = glob(f"{self.folder}/*_p?_f1*.log")
        assert (
            full_molecule_a
        ), f"No files named *_p?_f1.log in folder {self.folder}."
        full_molecule_a1 = full_molecule_a[0].split("/")[-1]
        return full_molecule_a1.split("_p")[0]

    def _get_all_files_full_molecule(self):
        all_files_full_molecule = []
        all_files_full_molecule_indices = []
        # find all files for full molecule
        full_molecule_pattern = re.compile(
            gaussian_dias_filename_point_without_fragment
        )

        # file all the files that match
        for filename in glob(f"{self.folder}/*.log"):
            match = full_molecule_pattern.match(filename)
            if match:
                number_after_p = int(match.group(1))
                all_files_full_molecule.append(filename)
                all_files_full_molecule_indices.append(number_after_p)

        # Sort filenames based on indices
        return sorted(
            all_files_full_molecule,
            key=lambda x: all_files_full_molecule_indices[
                all_files_full_molecule.index(x)
            ],
        )

    def _get_all_files_fragment1(self):
        all_files_fragment1 = []
        all_files_fragment1_indices = []

        fragment1_pattern = re.compile(
            gaussian_dias_filename_point_with_fragment1
        )

        for filename in glob(f"{self.folder}/*.log"):
            match = fragment1_pattern.match(filename)
            if match:
                number_after_p = int(match.group(1))
                all_files_fragment1.append(filename)
                all_files_fragment1_indices.append(number_after_p)

        # Sort filenames based on indices
        return sorted(
            all_files_fragment1,
            key=lambda x: all_files_fragment1_indices[
                all_files_fragment1.index(x)
            ],
        )

    def _get_all_files_fragment2(self):
        all_files_fragment2 = []
        all_files_fragment2_indices = []

        fragment2_pattern = re.compile(
            gaussian_dias_filename_point_with_fragment2
        )

        for filename in glob(f"{self.folder}/*.log"):
            match = fragment2_pattern.match(filename)
            if match:
                number_after_p = int(match.group(1))
                all_files_fragment2.append(filename)
                all_files_fragment2_indices.append(number_after_p)

        # Sort filenames based on indices
        return sorted(
            all_files_fragment2,
            key=lambda x: all_files_fragment2_indices[
                all_files_fragment2.index(x)
            ],
        )

    def _get_all_files_reactants(self):
        all_files_reactants = []
        all_files_reactants_indices = []

        reactants_pattern = re.compile(gaussian_dias_filename_with_reactant)

        for filename in glob(f"{self.folder}/*.log"):
            match = reactants_pattern.match(filename)
            if match:
                number_after_r = int(match.group(1))
                all_files_reactants.append(filename)
                all_files_reactants_indices.append(number_after_r)

        # Sort filenames based on indices
        return sorted(
            all_files_reactants,
            key=lambda x: all_files_reactants_indices[
                all_files_reactants.index(x)
            ],
        )

    def _get_all_energies(self, list_of_files):
        all_energies = []
        for file in list_of_files:
            gout = Gaussian16Output(filename=file)
            energy = gout.molecule.energy
            # convert energy from eV (ase units) to kcal/mol
            energy /= units.kcal / units.mol
            all_energies.append(energy)
        return all_energies


class ORCADIASOutFolder(DIASOutputFolder):
    def __init__(
        self,
        folder,
        atom1,
        atom2,
        zero=False,
        ref_file=None,
        outputname="dias",
    ):
        super().__init__(
            folder=folder,
            atom1=atom1,
            atom2=atom2,
            zero=zero,
            ref_file=ref_file,
            outputname=outputname,
        )

    @property
    def job_basename(self):
        full_molecule_a = glob(f"{self.folder}/*_p?_f1*out")
        assert (
            full_molecule_a
        ), f"No files named *_p?_f1*out in folder {self.folder}."
        full_molecule_a1 = full_molecule_a[0].split("/")[-1]
        return full_molecule_a1.split("_p")[0]

    def _get_all_files_full_molecule(self):
        all_files_full_molecule = []
        all_files_full_molecule_indices = []
        # find all files for full molecule
        full_molecule_pattern = re.compile(
            orca_dias_filename_point_without_fragment
        )

        # file all the files that match
        for filename in glob(f"{self.folder}/*.out"):
            match = full_molecule_pattern.match(filename)
            if match:
                number_after_p = int(match.group(1))
                all_files_full_molecule.append(filename)
                all_files_full_molecule_indices.append(number_after_p)

        # Sort filenames based on indices
        return sorted(
            all_files_full_molecule,
            key=lambda x: all_files_full_molecule_indices[
                all_files_full_molecule.index(x)
            ],
        )

    def _get_all_files_fragment1(self):
        all_files_fragment1 = []
        all_files_fragment1_indices = []

        fragment1_pattern = re.compile(orca_dias_filename_point_with_fragment1)

        for filename in glob(f"{self.folder}/*.out"):
            match = fragment1_pattern.match(filename)
            if match:
                number_after_p = int(match.group(1))
                all_files_fragment1.append(filename)
                all_files_fragment1_indices.append(number_after_p)

        # Sort filenames based on indices
        return sorted(
            all_files_fragment1,
            key=lambda x: all_files_fragment1_indices[
                all_files_fragment1.index(x)
            ],
        )

    def _get_all_files_fragment2(self):
        all_files_fragment2 = []
        all_files_fragment2_indices = []

        fragment2_pattern = re.compile(orca_dias_filename_point_with_fragment2)

        for filename in glob(f"{self.folder}/*.out"):
            match = fragment2_pattern.match(filename)
            if match:
                number_after_p = int(match.group(1))
                all_files_fragment2.append(filename)
                all_files_fragment2_indices.append(number_after_p)

        # Sort filenames based on indices
        return sorted(
            all_files_fragment2,
            key=lambda x: all_files_fragment2_indices[
                all_files_fragment2.index(x)
            ],
        )

    def _get_all_files_reactants(self):
        all_files_reactants = []
        all_files_reactants_indices = []

        reactants_pattern = re.compile(orca_dias_filename_with_reactant)

        for filename in glob(f"{self.folder}/*.out"):
            match = reactants_pattern.match(filename)
            if match:
                number_after_r = int(match.group(1))
                all_files_reactants.append(filename)
                all_files_reactants_indices.append(number_after_r)

        # Sort filenames based on indices
        return sorted(
            all_files_reactants,
            key=lambda x: all_files_reactants_indices[
                all_files_reactants.index(x)
            ],
        )

    def _get_all_energies(self, list_of_files):
        all_energies = []
        for file in list_of_files:
            oout = ORCAOutput(filename=file)
            energy = oout.molecule.final_energy
            # convert energy from eV (ase units) to kcal/mol
            energy /= units.kcal / units.mol
            all_energies.append(energy)
        logger.info(f"All energies: {all_energies}")
        return all_energies
