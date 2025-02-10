import os
from functools import cached_property


class XTBInput:
    def __init__(self, filename):
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File {filename} not found.")
        self.filename = filename

    @cached_property
    def contents(self):
        with open(self.filename) as f:
            return [line.strip() for line in f.readlines()]

    @cached_property
    def content_groups(self):
        groups = {}
        current_group = None
        for line in self.contents:
            if line.startswith("$"):
                parts = line.split(maxsplit=1)
                current_group = parts[0].strip("$")
                groups[current_group] = [parts[1]] if len(parts) > 1 else []
            elif current_group:
                groups[current_group].append(line)
        return groups

    @property
    def charge(self):
        if "chrg" in self.content_groups:
            return int(self.content_groups["chrg"][0])
        return 0

    @property
    def spin(self):
        if "spin" in self.content_groups:
            return float(self.content_groups["spin"][0])
        """Total spin from number of unpaired electrons"""
        """default Nalpha-Nbeta of the molecule is 0.0"""
        return 0.0

    @property
    def method(self):
        if "gfn" in self.content_groups:
            for line in self.content_groups["gfn"]:
                if "method" in line:
                    if int(line.split("=")[1].strip()) == 0:
                        return "GFN0-xTB"
                    elif int(line.split("=")[1].strip()) == 1:
                        return "GFN1-xTB"
                    elif int(line.split("=")[1].strip()) == 2:
                        return "GFN2-xTB"
        return "GFN2-xTB"

    @property
    def opt_level(self):
        if "opt" in self.content_groups:
            for line in self.content_groups["opt"]:
                if "optlevel" in line:
                    return line.split("=")[1].strip()
        return "normal"

    @property
    def potential(self):
        if "wall" in self.content_groups:
            for line in self.content_groups["wall"]:
                if "potential" in line:
                    return line.split("=")[1].strip()
        return "polynomial"

    @property
    def electronic_temperature(self):
        """electronic temperature for the Fermi smearing in K"""
        if "scc" in self.content_groups:
            for line in self.content_groups["scc"]:
                if "temp" in line:
                    return float(line.split("=")[1].strip())
        return 300.0

    @property
    def scc_maxiterations(self):
        if "scc" in self.content_groups:
            for line in self.content_groups["scc"]:
                if "iterations" in line:
                    return int(line.split("=")[1].strip())
        return 250

    @property
    def engine(self):
        if "opt" in self.content_groups:
            for line in self.content_groups["opt"]:
                if "engine" in line:
                    return line.split("=")[1].strip()
        return "rf"

    @property
    def solvent(self):
        if "gbsa" in self.content_groups:
            for line in self.content_groups["gbsa"]:
                if "solvent" in line:
                    return line.split("=")[1].strip()
        return None

    @property
    def thermo_temperature(self):
        """temperature for thermostatistical calculation in K"""
        if "thermo" in self.content_groups:
            for line in self.content_groups["thermo"]:
                if "temp" in line:
                    return float(line.split("=")[1].strip())
        return 300.0

    def _get_md_setting(self, key, default):
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if key in line:
                    try:
                        return float(line.split("=")[1].strip())
                    except ValueError:
                        return default
            return default
        return None

    @property
    def md_dump(self):
        """interval for trajectory printout in fs"""
        return self._get_md_setting("dump", 50.0)

    @property
    def md_hmass(self):
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if "hmass" in line:
                    """multiple of the mass of hydrogen atoms"""
                    return int(line.split("=")[1].strip())
            return 4
        return None

    @property
    def md_velo(self):
        """control whether to write out velocities"""
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if "velo" in line:
                    if line.split("=")[1].strip().lower() in {"1", "true"}:
                        return True
            return False
        return None

    @property
    def md_nvt(self):
        """check whether a simulation in NVT ensemble is performed"""
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if "nvt" in line:
                    if line.split("=")[1].strip().lower() in {"0", "false"}:
                        return False
            return True
        return None

    @property
    def md_temperature(self):
        """thermostat temperature in K"""
        return self._get_md_setting("temp", 298.15)

    @property
    def md_time(self):
        """total run time of simulation in ps"""
        return self._get_md_setting("time", 50.0)

    @property
    def md_step(self):
        """time step for propagation in fs"""
        return self._get_md_setting("step", 4.0)