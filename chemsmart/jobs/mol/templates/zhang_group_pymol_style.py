# Xinglong Zhang Group PyMOL Style Settings
# Add_VDW creates a copy of an object with full-sized, transparent spheres
# Bondi VDW values added below to override default PyMOL settings

from pymol import cmd

# Bondi VDW values
cmd.alter("elem Ac", "vdw=2.00")
cmd.alter("elem Al", "vdw=2.00")
cmd.alter("elem Am", "vdw=2.00")
cmd.alter("elem Sb", "vdw=2.00")
cmd.alter("elem Ar", "vdw=1.88")
cmd.alter("elem As", "vdw=1.85")
cmd.alter("elem At", "vdw=2.00")
cmd.alter("elem Ba", "vdw=2.00")
cmd.alter("elem Bk", "vdw=2.00")
cmd.alter("elem Be", "vdw=2.00")
cmd.alter("elem Bi", "vdw=2.00")
cmd.alter("elem Bh", "vdw=2.00")
cmd.alter("elem B ", "vdw=2.00")
cmd.alter("elem Br", "vdw=1.85")
cmd.alter("elem Cd", "vdw=1.58")
cmd.alter("elem Cs", "vdw=2.00")
cmd.alter("elem Ca", "vdw=2.00")
cmd.alter("elem Cf", "vdw=2.00")
cmd.alter("elem C ", "vdw=1.70")
cmd.alter("elem Ce", "vdw=2.00")
cmd.alter("elem Cl", "vdw=1.75")
cmd.alter("elem Cr", "vdw=2.00")
cmd.alter("elem Co", "vdw=2.00")
cmd.alter("elem Cu", "vdw=1.40")
cmd.alter("elem Cm", "vdw=2.00")
cmd.alter("elem Ds", "vdw=2.00")
cmd.alter("elem Db", "vdw=2.00")
cmd.alter("elem Dy", "vdw=2.00")
cmd.alter("elem Es", "vdw=2.00")
cmd.alter("elem Er", "vdw=2.00")
cmd.alter("elem Eu", "vdw=2.00")
cmd.alter("elem Fm", "vdw=2.00")
cmd.alter("elem F ", "vdw=1.47")
cmd.alter("elem Fr", "vdw=2.00")
cmd.alter("elem Gd", "vdw=2.00")
cmd.alter("elem Ga", "vdw=1.87")
cmd.alter("elem Ge", "vdw=2.00")
cmd.alter("elem Au", "vdw=1.66")
cmd.alter("elem Hf", "vdw=2.00")
cmd.alter("elem Hs", "vdw=2.00")
cmd.alter("elem He", "vdw=1.40")
cmd.alter("elem Ho", "vdw=2.00")
cmd.alter("elem In", "vdw=1.93")
cmd.alter("elem I ", "vdw=1.98")
cmd.alter("elem Ir", "vdw=2.00")
cmd.alter("elem Fe", "vdw=2.00")
cmd.alter("elem Kr", "vdw=2.02")
cmd.alter("elem La", "vdw=2.00")
cmd.alter("elem Lr", "vdw=2.00")
cmd.alter("elem Pb", "vdw=2.02")
cmd.alter("elem Li", "vdw=1.82")
cmd.alter("elem Lu", "vdw=2.00")
cmd.alter("elem Mg", "vdw=1.73")
cmd.alter("elem Mn", "vdw=2.00")
cmd.alter("elem Mt", "vdw=2.00")
cmd.alter("elem Md", "vdw=2.00")
cmd.alter("elem Hg", "vdw=1.55")
cmd.alter("elem Mo", "vdw=2.00")
cmd.alter("elem Nd", "vdw=2.00")
cmd.alter("elem Ne", "vdw=1.54")
cmd.alter("elem Np", "vdw=2.00")
cmd.alter("elem Ni", "vdw=1.63")
cmd.alter("elem Nb", "vdw=2.00")
cmd.alter("elem N ", "vdw=1.55")
cmd.alter("elem No", "vdw=2.00")
cmd.alter("elem Os", "vdw=2.00")
cmd.alter("elem O ", "vdw=1.52")
cmd.alter("elem Pd", "vdw=1.63")
cmd.alter("elem P ", "vdw=1.80")
cmd.alter("elem Pt", "vdw=1.72")
cmd.alter("elem Pu", "vdw=2.00")
cmd.alter("elem Po", "vdw=2.00")
cmd.alter("elem K ", "vdw=2.75")
cmd.alter("elem Pr", "vdw=2.00")
cmd.alter("elem Pm", "vdw=2.00")
cmd.alter("elem Pa", "vdw=2.00")
cmd.alter("elem Ra", "vdw=2.00")
cmd.alter("elem Rn", "vdw=2.00")
cmd.alter("elem Re", "vdw=2.00")
cmd.alter("elem Rh", "vdw=2.00")
cmd.alter("elem Rb", "vdw=2.00")
cmd.alter("elem Ru", "vdw=2.00")
cmd.alter("elem Rf", "vdw=2.00")
cmd.alter("elem Sm", "vdw=2.00")
cmd.alter("elem Sc", "vdw=2.00")
cmd.alter("elem Sg", "vdw=2.00")
cmd.alter("elem Se", "vdw=1.90")
cmd.alter("elem Si", "vdw=2.10")
cmd.alter("elem Ag", "vdw=1.72")
cmd.alter("elem Na", "vdw=2.27")
cmd.alter("elem Sr", "vdw=2.00")
cmd.alter("elem S ", "vdw=1.80")
cmd.alter("elem Ta", "vdw=2.00")
cmd.alter("elem Tc", "vdw=2.00")
cmd.alter("elem Te", "vdw=2.06")
cmd.alter("elem Tb", "vdw=2.00")
cmd.alter("elem Tl", "vdw=1.96")
cmd.alter("elem Th", "vdw=2.00")
cmd.alter("elem Tm", "vdw=2.00")
cmd.alter("elem Sn", "vdw=2.17")
cmd.alter("elem Ti", "vdw=2.00")
cmd.alter("elem W ", "vdw=2.00")
cmd.alter("elem U ", "vdw=1.86")
cmd.alter("elem V ", "vdw=2.00")
cmd.alter("elem Xe", "vdw=2.16")
cmd.alter("elem Yb", "vdw=2.00")
cmd.alter("elem Y ", "vdw=2.00")
cmd.alter("elem Zn", "vdw=1.39")
cmd.alter("elem Zr", "vdw=2.00")
cmd.rebuild()


def default_render():
    """Set default rendering settings for PyMOL."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.util.ray_shadows("light")
    cmd.set("specular", 0.25)
    cmd.set("spec_power", 300)
    cmd.set("spec_reflect", 0.5)
    cmd.set("antialias", 1)
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 45)
    cmd.set("transparency", 0)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.3)


cmd.extend("default_render", default_render)


def movie_render():
    """Set rendering settings for PyMOL movie."""
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")
    cmd.set("specular", 0.25)
    cmd.set("spec_power", 300)
    cmd.set("spec_reflect", 0.5)
    cmd.util.ray_shadows("none")

    cmd.set("antialias", 1)
    cmd.set("orthoscopic", 0)
    cmd.set("field_of_view", 45)
    cmd.set("transparency", 0)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_frame", 1)


cmd.extend("movie_render", movie_render)


def enhance_visuals():
    # Set background color to white
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", "off")

    # Enable depth cueing (fog effect)
    cmd.set("depth_cue", 1)
    cmd.set("fog_start", 0.3)

    # Adjust lighting and material properties
    cmd.set("ambient", 0.5)
    cmd.set("spec_count", 5)
    cmd.set("shininess", 50)
    cmd.set("specular", 1)
    cmd.set("reflect", 0.1)

    # Set antialiasing for smoother visuals
    cmd.set("antialias", 3)

    # Set ray tracing mode for high-quality rendering
    cmd.set("ray_trace_mode", 1)


# Extend PyMOL with the new function
cmd.extend("enhance_visuals", enhance_visuals)


def ballnstick(arg1):
    """Set ball-and-stick representation for a given selection."""
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray85", "elem C and " + arg1)
    cmd.color("gray98", "elem H and " + arg1)
    cmd.color("gray", "elem Ag and " + arg1)
    cmd.set("stick_radius", 0.07, arg1)
    cmd.set("sphere_scale", 0.15, arg1)
    cmd.alter(arg1 + " and elem H", "vdw=0.75")
    cmd.set("stick_color", "black", arg1)
    cmd.set("dash_gap", 0.2)
    cmd.set("dash_radius", 0.05)
    cmd.set("label_distance_digits", 2)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)


def pymol_style(arg1):
    default_render()
    ballnstick(arg1)


cmd.extend("pymol_style", pymol_style)


def cylview_style(arg1):
    enhance_visuals()
    ballnstick(arg1)


cmd.extend("cylview_style", cylview_style)


def movie_style(arg1):
    movie_render()
    ballnstick(arg1)


cmd.extend("movie_style", movie_style)


def add_vdw(arg1):
    cmd.copy(arg1 + "_vdw", arg1)
    cmd.alter(arg1 + "_vdw and elem H", "vdw=1.09")
    cmd.rebuild()
    cmd.set("sphere_scale", 1, arg1 + "_vdw")
    cmd.hide("nonbonded", arg1 + "_vdw")
    cmd.hide("lines", arg1 + "_vdw")
    cmd.hide("sticks", arg1 + "_vdw")
    cmd.set("sphere_transparency", 0.7, arg1 + "_vdw")


cmd.extend("add_vdw", add_vdw)


def nci(arg1, isosurface=0.5, range=1.0):
    dens_file = arg1 + "-dens"
    grad_file = arg1 + "-grad"
    cmd.isosurface("grad", grad_file, isosurface)
    cmd.ramp_new("ramp", dens_file, [-range, 0, range], "rainbow")
    cmd.set("surface_color", "ramp", "grad")
    cmd.set("two_sided_lighting", value=1)
    cmd.set("transparency", 0.5)
    cmd.set("surface_quality", 1)


cmd.extend("nci", nci)


def nci_binary(arg1, isosurface=0.5, range=1.0):
    dens_file = arg1 + "-dens"
    grad_file = arg1 + "-grad"
    cmd.isosurface("grad", grad_file, isosurface)
    cmd.ramp_new("ramp", dens_file, [-range, 0, range], "[blue,white,red]")
    cmd.set("surface_color", "ramp", "grad")
    cmd.set("two_sided_lighting", value=1)
    cmd.set("transparency", 0.5)
    cmd.set("surface_quality", 1)


cmd.extend("nci_binary", nci_binary)
