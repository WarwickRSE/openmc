import matplotlib.pyplot as plt
import openmc
#Adapted from OpenMC docs

openmc.Materials.cross_sections = '/media/raid/MathRadData/endfb-viii.1-hdf5/cross_sections.xml'

# Create materials
materials = openmc.Materials()

water = openmc.Material()
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.set_density('g/cm3', 1.0)
water.temperature = 300 # K

materials.append(water)
materials.export_to_xml()

# Create geometry
geometry = openmc.Geometry()
box = openmc.model.RectangularPrism(width=1.0, height=1.0,
                                    boundary_type='reflective')
ghost = openmc.Cell(name='ghost')
ghost.fill = water
ghost.region = -box

universe = openmc.Universe(cells=[ghost])
geometry.root_universe = universe

geometry.export_to_xml()

geometry_plot = universe.plot()
fig = geometry_plot.figure
fig.savefig('tmp.png')


# Assign simulation settings
settings = openmc.Settings()

point = openmc.stats.Point((0, 0, 0))
source = openmc.IndependentSource(space=point)

settings.run_mode = 'fixed source'
settings.source = source
settings.batches = 100
settings.inactive = 10
settings.particles = 1000

settings.track = [
    (1, 1, 1),
    (1, 1, 2),
    (1, 1, 3),
]

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