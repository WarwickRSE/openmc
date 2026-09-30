import matplotlib.pyplot as plt
import openmc
import numpy as np

def fwhm_c(x, y, peak):
    hm = peak/2.0
    indexes = np.where(y > hm)[0]
    print(indexes)
    return x[indexes[-1]] - x[indexes[0]]


def fwhm_y_for_each_x(x_centers, y_centers, heating_y):
    """Return the FWHM along y for each x position in a 2D heating map."""
    fwhm = np.full(x_centers.shape[0], np.nan)
    for i, profile in enumerate(heating_y):
        if np.all(profile == 0.0):
            continue
        peak = profile.max()
        if peak <= 0.0:
            continue
        half_max = peak / 2.0
        idx = np.where(profile >= half_max)[0]
        if idx.size < 2:
            continue
        fwhm[i] = y_centers[idx[-1]] - y_centers[idx[0]]
    return fwhm


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
    extent=(0.0, 20.0, -1.0, 1.0),
    aspect="auto",
    cmap="inferno",
)

ax.set_xlabel("x [cm]")
ax.set_ylabel("y [cm]")
fig.colorbar(image, ax=ax, label="Heating [eV/source particle]")
fig.tight_layout()
fig.savefig("heating_xy.png", dpi=200)




x_centers = np.linspace(0.0, 20.0, mesh.dimension[0], endpoint=False)
x_centers += 0.5 * (20.0 / mesh.dimension[0])
heating_lineout = heating_data.sum(axis=(1, 2))

peak = heating_lineout.max()
peak_ind = heating_lineout.argmax(axis=0)
peak_x = x_centers[peak_ind]
fwhm = fwhm_c(x_centers, heating_lineout, peak)
print(peak_x, fwhm)

# Compute the heating profile in y for each x position and its FWHM.
# This keeps the z dimension summed so each x slice is a 1D profile in y.
y_profile = heating_data.sum(axis=2)
y_lower = mesh.lower_left[1]
y_upper = mesh.upper_right[1]
y_step = (y_upper - y_lower) / mesh.dimension[1]
y_centers = np.linspace(y_lower, y_upper, mesh.dimension[1], endpoint=False)
y_centers += 0.5 * y_step
fwhm_y = fwhm_y_for_each_x(x_centers, y_centers, y_profile)
print("x positions with finite y-FWHM:", np.isfinite(fwhm_y).sum())

fig, ax = plt.subplots(figsize=(12, 4))
finite = np.isfinite(fwhm_y)
ax.plot(x_centers[finite], fwhm_y[finite], color="tab:blue", marker="o", markersize=3, label="FWHM in y")

# Fit a smooth, non-linear trend to the measured FWHM values.
coeffs = np.polyfit(x_centers[finite], fwhm_y[finite], 3)
fit_y = np.polyval(coeffs, x_centers)
print(coeffs)
ax.plot(x_centers, fit_y, color="tab:red", linewidth=2, label="Cubic best-fit")

ax.set_xlabel("x [cm]")
ax.set_ylabel("FWHM in y [cm]")
ax.set_title("Heating FWHM in y for each x position")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig("heating_fwhm_y_vs_x.png", dpi=200)

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

