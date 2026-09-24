
#We're going to disable the angular effects, and test penetration depth into a material
# This means doing a moderate amount of protons and then finding the peak of the depth. This we check against expectation 

# Then we repeat with the straggling enabled and check the width

# We do these for water and for a block of e.g. carbon


import openmc
from pytest import approx
import numpy as np

def test_proton_penetration(run_in_tmpdir):
    # Data path
    proton_path = '/media/raid/MathRadData/protons/'
    openmc.config['proton_data'] = proton_path

    # All-in-one model, rather than XML
    model = openmc.Model()

    # Create materials
    materials = openmc.Materials()

    water = openmc.Material()
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    #water.add_nuclide('H2', 2.0)
    water.set_density('g/cm3', 1.0)
    #water.set_density('g/cm3', 0.0006) # Steam
    water.temperature = 300 # K
    water.mean_excitation_energy = openmc.proton_data.mean_excitation_energy("water", proton_path + "mean_excitation_energies.txt") #eV

    materials.append(water)

    xlen = 20.0

    # Create geometry
    xmin = openmc.XPlane(x0=0.0, boundary_type='vacuum')
    xvoid = openmc.XPlane(x0=1.0)
    xmax = openmc.XPlane(x0=xlen, boundary_type='vacuum')
    ymin = openmc.YPlane(y0=-1.0, boundary_type='vacuum')
    ymax = openmc.YPlane(y0=1.0, boundary_type='vacuum')
    zmin = openmc.ZPlane(z0=-1.0, boundary_type='vacuum')
    zmax = openmc.ZPlane(z0=1.0, boundary_type='vacuum')

    transverse_region = +ymin & -ymax & +zmin & -zmax
    void_cell = openmc.Cell(region=+xmin & -xvoid & transverse_region)
    water_cell = openmc.Cell(
        fill=water,
        region=+xvoid & -xmax & transverse_region,
    )
    model.geometry = openmc.Geometry([void_cell, water_cell])

    model.settings.source = openmc.IndependentSource(
        particle='proton',
        space=openmc.stats.Point((0.1, 0.0, 0.0)),
        angle=openmc.stats.Monodirectional((1.0, 0.0, 0.0)),
        energy=openmc.stats.Discrete([100.0e6], [1.0])
    ) 

    universe = openmc.Universe(cells=[void_cell, water_cell])
    model.geometry.root_universe = universe

    model.settings.run_mode = 'fixed source'
    model.settings.proton_transport = True # The default but showing that it is wired in
    model.settings.batches = 1
    model.settings.particles = 10

    model.settings.cutoff = {'energy_proton': 2.0e5} #Cutoff energy in eV
    model.settings.proton_settings = {'max_step_len': 0.2, 'min_step_len':0.05, 'max_energy_loss':1e6}
    # Similarly, these default to True, but are shown here for clarity
    model.settings.proton_settings['use_sph'] = False
    model.settings.proton_settings['use_large_angle'] = False
    model.settings.proton_settings['use_straggling'] = False

    mesh = openmc.RegularMesh()
    mesh.lower_left = (0.0, -1.0, -1.0)
    mesh.upper_right = (xlen, 1.0, 1.0)
    mesh.dimension = (500, 40, 40)

    heating = openmc.Tally(name="proton heating")
    heating.filters = [openmc.MeshFilter(mesh)]
    heating.scores = ["heating"]

    model.tallies = openmc.Tallies([heating])

    model.run(apply_tally_results=True)
    
    print(heating)
    heating_data = heating.get_reshaped_data(expand_dims=True).squeeze()
    
    heating_lineout = heating_data.sum(axis=(1, 2))
    peak_ind = heating_lineout.argmax(axis=0)
    x_centers = np.linspace(0.0, xlen, mesh.dimension[0], endpoint=False)
    x_centers += 0.5 * (xlen / mesh.dimension[0])
    peak_x = x_centers[peak_ind]
    
    assert peak_x == approx(7.78 + 1)  # 1cm vacuum block at start, so 1cm 
