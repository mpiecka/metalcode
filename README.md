# metalcode
Automatic tool focused on deriving metallicities of open clusters. Based on the method described in Pöhnl &amp; Paunzen (2010, https://ui.adsabs.harvard.edu/abs/2010A%26A...514A..81P/abstract).


## Description
This is the version 1.0 of the automated version of the procedure devised by Pöhnl &amp; Paunzen (2010). The tool is focused on calculating metallicities Z (and logAge) of open clusters, assuming that accurate values of reddening and distance are provided (determined by independent methods). Before the code is applied, data for the cluster members (photometric brightness and colour) need to be prepared together with a file containing the list of clusters together with additional parameters. Examples of the data files are provided with the code. The code is applicable to Johnson (V,B-V), Gaia (G,BP-RP) and 2MASS (J,J-Ks) photometric systems.

The run-time of the code will depend on the total number of clusters in the included list, on the number of cluster members, and on the user input parameters. For example, if user inputs:

* delta_logAge = 0.2
* delta_Z = 0.005
* Nredd = 5
* Niter = 6

code will return results within 1-2 min for a typical open cluster. However, it can run for longer in the case of a larger cluster (for the included example of NGC 6791, the code returned results after 20 min). Furthermore, the specific run-time is also hardware-dependent (the code was tested on AMD Ryzen 3 PRO 4450U).

See Piecka & Paunzen (submitted) for a full description of the methods applied in the code.


## Requirements
In order to run the code, user must have installed Python 3 with `numpy`. The other libraries (`matplotlib`, `time`, `os`) are not required for the proper functionality of the code, but provide additional information useful information (e.g. figures).

The code was tested on the following operating systems:

* Windows 10
* Ubuntu 20.04 LTS
* Fedora 34


## Installation
Only Python 3 and the mentioned libraries need to be installed. Otherwise, no additional installation is needed.


## Usage
To launch the tool, run the script `metalcode_v1_0.py`. For successful application of the tool, a cluster list and the associated data files need to be included prior to running the script.


## Input
We describe several data files in this section of the documentation. As column separation, we use spaces between values. Furthermore, isochrone grids are required for the code to run. The included grids (logAge=6.6..10.0, Z=0.005..0.040, delta_logAge=0.1, delta_Z=0.005) are for the three photometric systems described below. The isochrones should be included in the main folder, the other files (described below) should be located in the `clusters` folder.

On the input, user must provide a file containing the list of clusters together with additional parameters. The structure of the file adheres to the following format (the first line of the file is skipped on loading):

```
CLUSTER_NAME   GAL_LATITUDE_deg   PARALLAX_mas   DISTANCE_pc   E(B-V)_mag
...            ...                ...            ...           ...
```

The cluster name should be written as one word (spaces should be replaced by underscores). Galactic latitude and parallax are not necessary - they should be used only if reddening is taken from extinction maps (in that case, `expcor` parameter in the code should be changed to 1). If the reddening value is not known and there is no good guess, set the value to be any negative value. The code will then use a pre-determined set of reddening values (in magnitudes: 0.010, 0.040, 0.080, 0.125, 0.250, 0.500, 0.750, 1.000, 1.500, 2.000).

Secondly, a set of files containing cluster data is required. The cluster data should be provided for the specific photometric system, and the file name should coincide with `CLUSTER_NAME_X`, where the suffix `X` should be replaced by the following:

* `G` for Gaia (G, BP-RP)
* `2` for 2MASS (J, J-Ks)
* `J` for Johnson (V, B-V)
 
The first line of the data file is skipped. The columns should follow the given format (in mag):

```
PHOTOMETRIC_BRIGHTNESS   PHOTOMETRIC_COLOUR
...                      ...
```

We strongly suggest that the users pre-analyse the colour-magnitude diagrams. Obvious binary sequences, white dwarfs, and possible other clear outliers should be removed in advance. This is necessary in the current version of the code due to the limitations of the included isochrone fitting sub-procedure.

Finally, the code will ask the user to specify additional parameters once it has
been launched.

1. **Photometric system:** Enter G, J or 2 (depending on the photometric system
   for which the data are available; see Photometric System table).
2. **Isochrone grid spacing, age:** In the current version, the user can choose
   between two spacings in the isochrone grid (0.1 or 0.2).
3. **Isochrone grid spacing, Z:** In the current version, use only values 0.005
   (can be changed by the user, but the set of isochrones should be changed
   accordingly, if necessary).
4. **Number of reddening iterations:** The number of reddening values that
   should be studied by the code. Choose 1 if you want to study only the initial
   estimate value. For 0, a predetermined set of ten values is used. Otherwise,
   use an odd number.
5. **Reddening range:** The relative range for reddening iterations. For
   example, if 0.3 is given and `Nredd` > 1, then the code will start at the
   value `0.7*E(B-V)_ini` and ends at `1.3*E(B-V)`. The value of the initial
   estimate is always included (if `Nredd` >= 1). Values between 0 and 1 are
   acceptable, excluding the limits.
6. **Maximum number of iterations:** Determines the maximum number of iterations
   while searching for metallicity for a given reddening value. Necessary
   because the code may get stuck between two possible solutions. Should not be
   a large number, because the number of iterations is typically smaller than
   five. We recommend using 6 for the currently included grids.
