def CombineH5Files(path):
    import os
    import h5py
    import numpy as np

    steps = [int(p.split(".")[1]) for p in os.listdir(path) if "params" in p]
    steps.sort()
    goodsteps = []

    with h5py.File(f"{path}/data.h5", "w") as file:
        for step in steps:
            # read/write attrs
            with h5py.File(f"{path}/params.{step:05d}", "r") as params:
                if step == steps[0]:
                    for attr in params.keys():
                        file.attrs[attr] = params[attr][()][0]
            # read/write specs
            try:
                with h5py.File(f"{path}/spec/spec.tot.{step:05d}", "r") as specfile:
                    spectrum = {
                        k: np.sum(specfile[k][:], axis=(1, 2, 3))
                        if not k.startswith("nr")
                        else specfile[k][:]
                        for k in specfile.keys()
                        if k.startswith("n")
                    }
                    spectrum["e"] = specfile["ebins"][:]
                    for k in spectrum.keys():
                        file.create_dataset(
                            f"Step{step:05d}/spectra/{k}",
                            data=spectrum[k],
                            chunks=True,
                        )
                    file[f"Step{step:05d}"].attrs["time"] = (
                        file.attrs["output:start"]
                        + file.attrs["output:interval"] * step
                    )
                goodsteps.append(step)
            except Exception as e:
                print("Unable to open file for step", step, e)
        file.attrs["NumSteps"] = len(goodsteps)


def get(fname):
    """
    Reads in a hdf5 file at `fname`, extracts the spectra and metadata, and returns them in a lazy xarray dataset.
    Parameters
    ----------
    fname : str
        The file name of the hdf5 file.
    Returns
    -------
    xarray.Dataset
        A lazy-loaded xarray dataset containing the spectra and metadata read from the hdf5 file.
    """
    import h5py
    import dask.array as da
    import numpy as np
    import xarray as xr

    file = h5py.File(fname, "r")
    step0 = list(file.keys())[0]
    nsteps = file.attrs["NumSteps"]
    times = np.array([file[f"Step{s:05d}"].attrs["time"][()] for s in range(nsteps)])

    ds = xr.Dataset()

    spectra = [s for s in file[step0]["spectra"].keys() if s.startswith("n")]

    for k in file.attrs.keys():
        if type(file.attrs[k]) == bytes or type(file.attrs[k]) == np.bytes_:
            ds.attrs[k] = file.attrs[k].decode("UTF-8")
        else:
            ds.attrs[k] = file.attrs[k]

    for s in spectra:
        dask_arrays = []
        for st in range(nsteps):
            array = da.from_array(file[f"Step{st:05d}/spectra/{s}"])
            dask_arrays.append(array[:])
        x = xr.DataArray(
            da.stack(dask_arrays, axis=0),
            dims=["t", "e"],
            name=s,
            coords={
                "t": times,
                "e": file[f"Step{0:05d}/spectra/e"][:],
            },
        )
        ds[s] = x
    return ds
