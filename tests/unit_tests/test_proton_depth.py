
#We're going to disable the angular effects, and test penetration depth into a material
# This means doing a moderate amount of protons and then finding the peak of the depth. This we check against expectation 

# Then we repeat with the straggling enabled and check the width

# We do these for water and for a block of e.g. carbon


import openmc
import pytest
import numpy as np
from os import environ

def create_simple_water():
    """ Water with standard element mix """
    proton_path = environ['OPENMC_PROTON_DATA']
     # Create materials

    water = openmc.Material()
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.set_density('g/cm3', 1.0)
    water.temperature = 300 # K
    water.mean_excitation_energy = openmc.proton_data.mean_excitation_energy("water", proton_path + "/mean_excitation_energies.txt") #eV

    return water

def create_simple_carbon():
    """ Pure carbon at 2g/cc """
    proton_path = environ['OPENMC_PROTON_DATA']
     # Create materials

    carbon = openmc.Material()
    carbon.add_element('C', 1.0)
    carbon.set_density('g/cm3', 2.0)
    carbon.temperature = 300 # K
    carbon.mean_excitation_energy = openmc.proton_data.mean_excitation_energy("carbon", proton_path + "/mean_excitation_energies.txt") #eV

    return carbon

def create_bone():
    """Material representing bone, matched to SDE PoC code """
    
    proton_path = environ['OPENMC_PROTON_DATA']

    bone = openmc.Material()
    bone.add_nuclide('H1', 0.0654709, percent_type='wo')
    bone.add_nuclide('C12', 0.536944, percent_type='wo')
    bone.add_nuclide('N14', 0.0215, percent_type='wo')
    bone.add_nuclide('O16', 0.032085, percent_type='wo')
    bone.add_nuclide('F19', 0.167411, percent_type='wo')
    bone.add_nuclide('Ca40', 0.176589, percent_type='wo')

    bone.set_density('g/cm3', 1.45)
    bone.mean_excitation_energy = openmc.proton_data.mean_excitation_energy("bone", proton_path + "/mean_excitation_energies.txt") #eV
    return bone

def create_block(material, length):
    """Create a block of specified material, [0,length] x [-1,1] x [-1,1] width""" 

    materials = openmc.Materials()

    if(material == 'water'):
        fill = create_simple_water()
        materials.append(fill)
    elif(material == "carbon"):
        fill = create_simple_carbon()
        materials.append(fill)
    elif(material == "bone"):
        fill = create_bone()
        materials.append(fill)
    else:
        raise ValueError("Material not known")

    # Create geometry
    xmin = openmc.XPlane(x0=0.0, boundary_type='vacuum')
    xmax = openmc.XPlane(x0=length, boundary_type='vacuum')
    ymin = openmc.YPlane(y0=-1.0, boundary_type='vacuum')
    ymax = openmc.YPlane(y0=1.0, boundary_type='vacuum')
    zmin = openmc.ZPlane(z0=-1.0, boundary_type='vacuum')
    zmax = openmc.ZPlane(z0=1.0, boundary_type='vacuum')

    transverse_region = +ymin & -ymax & +zmin & -zmax
    mat_cell = openmc.Cell(
        fill=fill,
        region=+xmin & -xmax & transverse_region,
    )
    geometry = openmc.Geometry([mat_cell])

    universe = openmc.Universe(cells=[mat_cell])
    geometry.root_universe = universe

    return(materials, geometry)

def water_block():
    """ Create a block of SDP water of length 20cm"""
    model = openmc.Model()
    model.materials, model.geometry = create_block('water', 20.0)
    return model

def carbon_block():
    """Create a block of 2g/cc carbon of length 20cm"""
    model = openmc.Model()
    model.materials, model.geometry = create_block('carbon', 20.0)
    return model

def bone_block():
    """Create a block compact bone of length 30cm"""
    model = openmc.Model()
    model.materials, model.geometry = create_block('bone', 30.0)
    return model

def run_model(model, xlen, energy, expected_peak):

    model.settings.run_mode = 'fixed source'
    model.settings.proton_transport = True # The default but showing that it is wired in
    model.settings.batches = 1
    model.settings.particles = 10
    model.settings.verbosity = 1

    source = openmc.IndependentSource(
        particle='proton',
        space=openmc.stats.Point((0.001, 0.0, 0.0)),
        angle=openmc.stats.Monodirectional((1.0, 0.0, 0.0)),
        energy=openmc.stats.Discrete([energy], [1.0])
    )

    model.settings.source = source
    model.settings.cutoff = {'energy_proton': 2.0e5} #Cutoff energy in eV
    model.settings.proton_settings = {'max_step_len': 0.2, 'min_step_len':0.05, 'max_energy_loss':1e6}
    # Similarly, these default to True, but are shown here for clarity
    model.settings.proton_settings['use_sph'] = False
    model.settings.proton_settings['use_large_angle'] = False
    model.settings.proton_settings['use_straggling'] = False

    mesh = openmc.RegularMesh()
    mesh.lower_left = (0.0, -1.0, -1.0)
    mesh.upper_right = (xlen, 1.0, 1.0)
    mesh.dimension = (500, 2, 2)

    heating = openmc.Tally(name="proton heating")
    heating.filters = [openmc.MeshFilter(mesh)]
    heating.scores = ["heating"]

    model.tallies = openmc.Tallies([heating])

    model.run(apply_tally_results=True)
    
    heating_data = heating.get_reshaped_data(expand_dims=True).squeeze()
    heating_lineout = heating_data.sum(axis=(1, 2))
    peak_ind = heating_lineout.argmax(axis=0)
    x_centers = np.linspace(0.0, xlen, mesh.dimension[0], endpoint=False)
    x_centers += 0.5 * (xlen / mesh.dimension[0])
    peak_x = x_centers[peak_ind]
    
    assert peak_x == pytest.approx(expected_peak, 1e-2)

def test_proton_penetration_water(run_in_tmpdir):
    """Penetration into water at SDP for 100MeV"""
    model = water_block()
    xlen = 20.0  # cm
    expected_peak = 7.72  #cm https://physics.nist.gov/cgi-bin/Star/ap_table.pl
    energy = 100e6 # MeV
    run_model(model, xlen, energy, expected_peak)

def test_proton_penetration_carbon(run_in_tmpdir):
    """Penetration into 2g/cc pure carbon at 150MeV"""
    model = carbon_block()
    xlen = 20.0
    expected_peak = 17.7  # I think this is right, from https://physics.nist.gov/cgi-bin/Star/ap_table.pl
    # TODO - compare to test code
    energy = 150e6
    run_model(model, xlen, energy, expected_peak)

@pytest.mark.parametrize('index', [50, 100, 150, 200])
def test_proton_penetration_bone(run_in_tmpdir, index):
    """Penetration into 'bone' block at 50, 100, 150 and 200 MeV"""
    model = bone_block()
    xlen = 30.0
     # I think this is right, from https://physics.nist.gov/cgi-bin/Star/ap_table.pl
    energy_depths = {50:2.41, 100:8.32, 150:17.0, 200:27.9}
    #for energy, depth in energy_depths.items():
    #    run_model(model, xlen, energy*1e6, depth)
    run_model(model, xlen, index*1e6, energy_depths[index])

