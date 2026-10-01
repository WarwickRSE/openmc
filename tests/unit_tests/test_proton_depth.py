# Testing the energy loss elements of the proton transport model
# We check overall penetration with straggling and all angular effects disabled
# and compare overall penetration to NIST data
# Then we check FWHM with straggling enabled and compare to ...???

# These are NOT fast tests as we need to transport and score a moderate number of protons

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

def create_artificial_water():
    """ Water with standard element mix """
    proton_path = environ['OPENMC_PROTON_DATA']
     # Create materials

    water = openmc.Material()
    water.add_element('H', 2.0)
    water.add_element('O', 1.0)
    water.set_density('g/cm3', 1.0)
    water.temperature = 300 # K
    water.mean_excitation_energy = 1e6 #eV - Extremely high

    return water

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
    elif(material == "fake_water"):
        fill = create_artificial_water()
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

def artificial_block():
    """Create a block of a fake high I material, length 20cm"""
    #This is like water in terms of density etc, but has I set high so that the energy loss is very slow
    model = openmc.Model()
    model.materials, model.geometry = create_block('fake_water', 20.0)
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

def fwhm_c(x, y, peak):
    hm = peak/2.0
    indexes = np.where(y > hm)[0]

    return x[indexes[-1]] - x[indexes[0]]

def expected_energy_variance_slope(material, energy_eV):
    """Return the fixed-energy straggling variance slope in eV^2/cm."""
    energy = energy_eV * 1.0e-6  # MeV
    proton_mass = 938.346  # MeV
    beta_sq = (2.0 * proton_mass + energy) * energy / (proton_mass + energy) ** 2
    gamma_sq = 1.0 / (1.0 - beta_sq)
    alpha = 1.0 / 137.0
    hbar = 4.136e-21 / (2.0 * np.pi)  # MeV s
    speed_of_light = 2.99792458e10  # cm/s
    prefactor = 4.0 * np.pi * (alpha * hbar * speed_of_light)**2
    energy_factor = gamma_sq * (1.0 - beta_sq / 2.0)

    electron_density = sum(
        atom_density * 1.0e24 * openmc.data.zam(name)[0]
        for name, atom_density in material.get_nuclide_atom_densities().items()
    )
    return prefactor * energy_factor * electron_density * 1.0e12

def run_model_straggled(model, xlen, energy, expected_slope):

    model.settings.run_mode = 'fixed source'
    model.settings.proton_transport = True # The default but showing that it is wired in
    model.settings.batches = 1
    model.settings.particles = 1000
    model.settings.track = [
        (1, 1, particle_id)
        for particle_id in range(1, model.settings.particles + 1)
    ]
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
    model.settings.proton_settings['use_straggling'] = True

    model.run()
    
    # Pool track energies into depth bins, then fit the variance over the
    # populated 1--xlen cm range, matching the analysis in water_block.py.
    n_depth_bins = 200
    depth_edges = np.linspace(0.0, xlen, n_depth_bins + 1)
    depth_centers = 0.5 * (depth_edges[:-1] + depth_edges[1:])
    energy_sum = np.zeros(n_depth_bins)
    energy_sum_sq = np.zeros(n_depth_bins)
    energy_count = np.zeros(n_depth_bins)

    tracks = openmc.Tracks('tracks.h5')
    assert len(tracks) == model.settings.particles
    for track in tracks:
        for particle_track in track.filter(particle='proton'):
            states = particle_track.states
            r = states['r']
            energies = states['E']
            segment_lengths = np.sqrt(
                np.diff(r['x'])**2 + np.diff(r['y'])**2 + np.diff(r['z'])**2
            )
            path_length = np.concatenate(([0.0], np.cumsum(segment_lengths)))
            bin_indices = np.digitize(path_length, depth_edges) - 1
            valid = (bin_indices >= 0) & (bin_indices < n_depth_bins)
            indices = bin_indices[valid]
            values = energies[valid]
            energy_sum += np.bincount(
                indices, weights=values, minlength=n_depth_bins
            )
            energy_sum_sq += np.bincount(
                indices, weights=values**2, minlength=n_depth_bins
            )
            energy_count += np.bincount(indices, minlength=n_depth_bins)

    has_data = energy_count > 0
    energy_mean = np.full(n_depth_bins, np.nan)
    energy_variance = np.full(n_depth_bins, np.nan)
    energy_mean[has_data] = energy_sum[has_data] / energy_count[has_data]
    energy_variance[has_data] = (
        energy_sum_sq[has_data] / energy_count[has_data]
        - energy_mean[has_data]**2
    )
    energy_variance[has_data] = np.maximum(energy_variance[has_data], 0.0)
    fit_mask = has_data & (depth_centers >= 0.0) & (depth_centers <= xlen)
    assert np.count_nonzero(fit_mask) >= 2
    fitted_slope = np.polyfit(
        depth_centers[fit_mask], energy_variance[fit_mask], 1
    )[0]
    assert fitted_slope == pytest.approx(expected_slope, rel=0.2), (
        f"Fitted variance slope {fitted_slope:.6e} eV^2/cm does not match "
        f"expected {expected_slope:.6e} eV^2/cm"
    )

@pytest.mark.parametrize('energy', [50, 100, 150, 200])
def test_proton_spread_highI(run_in_tmpdir, energy):
    """Penetration spread into a mock material with very high I, 100MeV."""
    model = artificial_block()
    xlen = 20.0  # cm
    energy = energy * 1e6
    expected_slope = expected_energy_variance_slope(
        model.materials[0], energy
    )
    run_model_straggled(model, xlen, energy, expected_slope)

#def test_proton_spread_highI2(run_in_tmpdir):
#    """Penetration spread into a mock material with very high I, 100MeV."""
#    model = artificial_block()

#    energy = 100e6 # MeV
#    xlen = 20.0  # cm
#    expected_slope = expected_energy_variance_slope(
#        model.materials[0], energy
#    )
#    run_model_straggled(model, xlen, energy, expected_slope)


#def test_proton_spread_water(run_in_tmpdir):
#    """Penetration spread into water at SDP for 100MeV"""
#    model = water_block()
#    xlen = 20.0  # cm
#    expected_peak = 7.72  #cm https://physics.nist.gov/cgi-bin/Star/ap_table.pl
#    energy = 100e6 # MeV
#    expected_fwhm = 0.823  # TODO this comes from the test code doing a first pass. 
#    run_model_straggled(model, xlen, energy, expected_peak, expected_fwhm)


