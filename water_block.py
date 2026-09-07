import matplotlib.pyplot as plt
import openmc
import numpy as np
#Adapted from OpenMC docs

openmc.Materials.cross_sections = '/media/raid/MathRadData/endfb-viii.1-hdf5/cross_sections.xml'

# Create materials
materials = openmc.Materials()

water = openmc.Material()
water.add_element('H', 2.0)
water.add_element('O', 1.0)
#water.add_nuclide('H2', 2.0)
water.set_density('g/cm3', 0.01)
#water.set_density('g/cm3', 0.0006) # Steam
water.temperature = 300 # K

materials.append(water)
materials.export_to_xml()

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
geometry = openmc.Geometry([void_cell, water_cell])

source = openmc.IndependentSource(
    particle='proton',
    space=openmc.stats.Point((0.1, 0.0, 0.0)),
    angle=openmc.stats.Monodirectional((1.0, 0.0, 0.0)),
    energy=openmc.stats.Discrete([100.0e6], [1.0])
) 

universe = openmc.Universe(cells=[void_cell, water_cell])
geometry.root_universe = universe

geometry.export_to_xml()

geometry_plot = universe.plot()
fig = geometry_plot.figure
fig.savefig('tmp.png')


# Assign simulation settings
settings = openmc.Settings()

settings.run_mode = 'fixed source'
settings.source = source
settings.batches = 1
settings.inactive = 10
settings.particles = 1000
settings.verbosity = 7

settings.track = [(1, 1, particle) for particle in range(1, int(settings.particles/10))]


mesh = openmc.RegularMesh()
mesh.lower_left = (0.0, -1.0, -1.0)
mesh.upper_right = (xlen, 1.0, 1.0)
mesh.dimension = (500, 40, 40)

heating = openmc.Tally(name="proton heating")
heating.filters = [openmc.MeshFilter(mesh)]
heating.scores = ["heating"]

tallies = openmc.Tallies([heating])
tallies.export_to_xml()

settings.export_to_xml()

openmc.run()
#openmc.run(mpi_args=['mpiexec', '-n', '4'])

tracks = openmc.Tracks('tracks.h5')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
for track in tracks:
    track.plot(ax)

ax.set_box_aspect((1, 1, 1))
fig.savefig('particle_tracks.png', dpi=200)

y_values = []
z_values = []

for track in tracks:
    for particle_track in track:
        states = particle_track.states
        y_values.extend(states["r"]["y"])
        z_values.extend(states["r"]["z"])

print("y range:", min(y_values), max(y_values))
print("z range:", min(z_values), max(z_values))
print("mesh dy:", 2.0 / mesh.dimension[1])
print("mesh dz:", 2.0 / mesh.dimension[2])

with openmc.StatePoint("statepoint.1.h5") as sp:
    tally = sp.get_tally(name="proton heating")

    print("estimator:", tally.estimator)
    print("mean shape:", tally.mean.shape)
    print("nonzero bins:", np.count_nonzero(tally.mean))

with openmc.StatePoint("statepoint.1.h5") as statepoint:
    tally = statepoint.get_tally(name="proton heating")
    heating_data = tally.get_reshaped_data(expand_dims=True).squeeze()

z_index = mesh.dimension[2] // 2
heatmap = heating_data.sum(axis=2)
#heatmap = heating_data[:, :, z_index]

fig, ax = plt.subplots(figsize=(12, 4))
image = ax.imshow(
    heatmap.T,
    origin="lower",
    extent=(0.0, xlen, -1.0, 1.0),
    aspect="auto",
    cmap="inferno",
)

ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
fig.colorbar(image, ax=ax, label="Heating [eV/source particle]")
fig.tight_layout()
fig.savefig("heating_xy.png", dpi=200)

x_centers = np.linspace(0.0, xlen, mesh.dimension[0], endpoint=False)
x_centers += 0.5 * (xlen / mesh.dimension[0])
heating_lineout = heating_data.sum(axis=(1, 2))

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(x_centers, heating_lineout, color="black")
ax.set_xlabel("x [cm]")
ax.set_ylabel("Heating integrated over y,z [eV/source particle]")
ax.set_title("Heating line-out summed over y and z")
#ax.set_yscale("log")
peak_heating = np.max(heating_lineout)
#ax.set_ylim(peak_heating / 10.0, peak_heating)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("heating_lineout_x.png", dpi=200)