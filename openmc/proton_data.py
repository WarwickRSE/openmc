from pathlib import Path

def mean_excitation_energy(
		name: str,
		filename: str) -> float:
	"""Return the mean excitation energy for a material or nuclide name.

	Parameters
	----------
	name : str
		Name in the first column of the data file.
	filename : str
		Path to the whitespace-delimited data file. 

	Returns
	-------
	float
		Value from the second column.

	Raises
	------
	KeyError
		If ``name`` is not present in the first column.
	"""
	path = Path(filename)
	with path.open(encoding="utf-8") as data_file:
		for line in data_file:
			fields = line.split()
			if not fields or fields[0].startswith("#"):
				continue
			if fields[0] == name:
				return float(fields[1])

	raise KeyError(f"No mean excitation energy found for {name!r}")
