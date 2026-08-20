import matplotlib.pyplot as plt
import openmc
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

# Create geometry
domain = openmc.model.RectangularParallelepiped(
    xmin=0.0, xmax=50.0,
    ymin=-1.0, ymax=1.0,
    zmin=-1.0, zmax=1.0,
    boundary_type='vacuum',
)

water_cell = openmc.Cell(fill=water, region=-domain)
geometry = openmc.Geometry([water_cell])

source = openmc.IndependentSource(
    particle='proton',
    space=openmc.stats.Point((0.1, 0.0, 0.0)),
    angle=openmc.stats.Monodirectional((1.0, 0.0, 0.0)),
    energy=openmc.stats.Discrete([1.0e6], [1.0])
) 

universe = openmc.Universe(cells=[water_cell])
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

settings.track = [(1, 1, particle) for particle in range(1, int(settings.particles/10))]

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