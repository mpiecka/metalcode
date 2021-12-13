# metalcode
Automatic tool focused on deriving metallicities of open clusters. Based on the method described in Pöhnl &amp; Paunzen (2010, https://ui.adsabs.harvard.edu/abs/2010A%26A...514A..81P/abstract).


## Description
This is the version 1.0 of the automated version of the procedure devised by Pöhnl &amp; Paunzen (2010). The tool is focused on calculating metallicities Z (and logAge) of open clusters, assuming that accurate values of reddening and distance are provided (determined by independent methods). Before the code is applied, data for the cluster members (photometric brightness and colour) need to be prepared together with a file containing the list of clusters together with additional parameters. Examples of the data files are provided with the code. The code is applicable to Johnson (V,B-V), Gaia (G,BP-RP) and 2MASS (J,J-Ks) photometric systems.

The run-time of the code will depend on the total number of clusters in the included list, on the number of cluster members, and on the user input parameters. For example, if user inputs:

* `age_step`=0.2
* `z_step`=0.005
* `Nredd`=5
* `Niter`=6

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

On the input (before the code is executed), the user must provide a file containing the list of clusters together with additional parameters (`\_complete.txt` in `clusters` folder). The structure of this file adheres to the following format (the first line of the file is skipped on loading):

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

1. **Photometric system**: Enter G, J or 2 (depending on the photometric system for which the data are available, see above for details).
2. **Isochrone grid spacing,** `age_step`: In the current version, the user can choose between two spacings in the isochrone grid (0.1 or 0.2).
3. **Isochrone grid spacing,** `z_step`: In the current version, use only value 0.005 (can be changed by the user, but the set of isochrones should be changed accordingly, if necessary).
4. **Number of reddening iterations**, `Nredd`: The number of reddening values that should be studied by the code. Choose 1 if you want to use only the initial estimate value `E(B-V)_ini`. For 0, a predetermined set of ten values is used. Otherwise, use any odd number larger than 1.
5. **Reddening range,** `redAdj`: The relative range for reddening iterations. For example, if `redAdj`=0.3 is given and `Nredd` > 1, then the code will start at the value `0.7*E(B-V)_ini` and end at `1.3*E(B-V)_ini`. The value of the initial estimate is always included (if `Nredd`>=1). Values between 0 and 1 are acceptable, excluding the limits.
6. **Maximum number of iterations,** `Niter`: Determines the maximum number of iterations while searching for metallicity for a given reddening value. Necessary because the code may get stuck between two possible solutions. A large number is not advised, because the number of iterations is typically smaller than five. We recommend using `Niter`=6 for the currently included grids.


## Output
The code provides all of the useful information on the output. If `debugTest` is set to True, the code will return additional information about the individual cluster members (values used in calculations, usually only required for debugging).

First of all, the solutions for different assumed reddening values will generally differ. For this, we include the results for all of the reddening values in a log-file in the `finished` folder. Included are the user input parameters, resulting cluster parameters (together with the quality-of-fit value, that should be minimised in the code) and the run-time for each of the individual clusters.

Secondly, the figures (CMD and LTN diagram) for the three best solutions are plotted saved in the `finished` folder. These figures should be consulted before interpreting the results.


## Sub-procedures
Details regarding the sub-procedures can be found in our paper. We would like to point out here that most of the sub-procedure can be easily exchanged. For example, the sub-procedures `metalcode_calib_absmg` and `metalcode_calib_clrex` are used to apply steps that deredden the colour and correct the brightness for the extinction. The transformation coefficients can be exchanged by the user (if required).

Furthermore, we use pre-prepared set of polynomial relation in order to calculate Teff and BC for a given combination of the colour and metallicity values. These calibrations were based on the isochrones themselves (and may slightly differ from the empirical, observation-based, relations found in the literature). If the user wishes to replace the relations, sets of polynomial coefficients have to be replaced in `metalcode_calib_tempe`. Because of how our code works, the user should prepare the coefficients for the different Z values, starting from Z=0.001 up to Z=0.040 (in the current version), with delta_Z=0.001.

Finally, the isochrone fitting technique is based only on a simple least-square method. In order to use any other technique, one should alter the file "metalcode_calc_lstsqr". The only requirement is that `LstSqr()` from this sub-procedure returns a quality-of-fit value that needs to be minimised.

We would like to point out that the currently included fitting technique was prepared only the for testing purposes, and it may not be sophisticated enough to produce results for proper scientific analysis. We urge the user to replace this sub-procedure if possible. In the future updates, we will replace this sub-procedure ourselves so that the code can be used for a scientific work right out of the box.


## Examples
We include a list of ten examples of open clusters that we analysed in our work. The observational data for the individual clusters were taken from the following sources:

* Gaia: cluster members (p>=0.70) from Cantat-Gaudin &amp; Anders (2020, https://ui.adsabs.harvard.edu/abs/2020A%26A...633A..99C/abstract). The same source was used to get the cluster parameters included in `\_complete.txt`, except for NGC 1039, where we had to use the reddening value from Dias et al. (2021, https://ui.adsabs.harvard.edu/abs/2021MNRAS.504..356D/abstract).
* 2MASS: the same cluster members as for Gaia. These stars were located using the positions (on the sky) from the list of the Gaia members.
* Johnson: we have used the data sets included in WEBDA.

All data files were manually pre-filtered in order to remove binary sequences, white dwarfs, and other possible outliers. A clear sequence of stars (main sequence + giants) is required with the currently introduced isochrone fitting sub-procedure.


## Acknowledgements
The work was supported from Operational Programme Research, Development and Education - ,,Project Internal Grant Agency of Masaryk University'' (No. CZ.02.2.69/0.0/0.0/19\_073/0016943).

This work makes use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.

This work makes use of data products from the Two Micron All Sky Survey, which is a joint project of the University of Massachusetts and the Infrared Processing and Analysis Center/California Institute of Technology, funded by the National Aeronautics and Space Administration and the National Science Foundation.

This research has made use of the WEBDA database (https://webda.physics.muni.cz), operated at the Department of Theoretical Physics and Astrophysics of the Masaryk University.

The isochrones were taken from http://stev.oapd.inaf.it/cgi-bin/cmd_3.5 (using default settings, except for the choice of the passbands).
