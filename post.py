import matplotlib.pyplot as plt
import openmc
import numpy as np




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

with openmc.StatePoint("statepoint.1.h5") as sp:
    tally = sp.get_tally(name="proton heating")

    mesh_filter = next(
        f for f in tally.filters
        if isinstance(f, openmc.MeshFilter)
    )
    mesh = mesh_filter.mesh

    print("mesh dimensions:", mesh.dimension)
    print("lower left:", mesh.lower_left)
    print("upper right:", mesh.upper_right)

    print("estimator:", tally.estimator)
    print("mean shape:", tally.mean.shape)
    print("nonzero bins:", np.count_nonzero(tally.mean))

    heating_data = tally.get_reshaped_data(expand_dims=True).squeeze()

z_index = mesh.dimension[2] // 2
heatmap = heating_data.sum(axis=2)
#heatmap = heating_data[:, :, z_index]

fig, ax = plt.subplots(figsize=(12, 4))
image = ax.imshow(
    heatmap.T,
    origin="lower",
    extent=(0.0, 50.0, -1.0, 1.0),
    aspect="auto",
    cmap="inferno",
)

ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
fig.colorbar(image, ax=ax, label="Heating [eV/source particle]")
fig.tight_layout()
fig.savefig("heating_xy.png", dpi=200)




x_centers = np.linspace(0.0, 60.0, mesh.dimension[0], endpoint=False)
x_centers += 0.5 * (60.0 / mesh.dimension[0])
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