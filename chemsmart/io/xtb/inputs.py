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
                if "method=" in line:
                    if int(line.split("=")[1]) == 1:
                        return "GFN1-xTB"
                    elif int(line.split("=")[1]) == 2:
                        return "GFN2-xTB"
        return "GFN2-xTB"

    @property
    def potential(self):
        if "wall" in self.content_groups:
            for line in self.content_groups["wall"]:
                if "potential=" in line:
                    return line.split("=")[1]
        return "polynomial"

    @property
    def engine(self):
        if "opt" in self.content_groups:
            for line in self.content_groups["opt"]:
                if "engine=" in line:
                    return line.split("=")[1]
        return "rf"

    @property
    def solvent(self):
        if "gbsa" in self.content_groups:
            for line in self.content_groups["gbsa"]:
                if "solvent=" in line:
                    return line.split("=")[1]
        return None

    def _get_md_setting(self, key, default):
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if f"{key}=" in line:
                    return float(line.split("=")[1])
            return default
        return None

    @property
    def md_dump(self):
        """interval for trajectory printout in fs"""
        return self._get_md_setting("dump", 50.0)

    @property
    def md_hydrogen_mass(self):
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if "hmass=" in line:
                    """multiple of the mass of hydrogen atoms"""
                    return int(line.split("=")[1])
            return 4
        return None

    @property
    def md_nvt(self):
        """check whether a simulation in NVT ensemble is performed"""
        if "md" in self.content_groups:
            for line in self.content_groups["md"]:
                if "nvt=" in line:
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