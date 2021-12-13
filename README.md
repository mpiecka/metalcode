# metalcode
Automatic tool focused on deriving metallicities of open clusters. Based on the method described in Pöhnl &amp; Paunzen (2010, https://ui.adsabs.harvard.edu/abs/2010A%26A...514A..81P/abstract).


## Description
This is the version 1.0 of the automated version of the procedure devised by Pöhnl &amp; Paunzen (2010). The tool is focused on calculating metallicities Z (and ages) of open clusters, assuming that accurate values of reddening and distance are provided (determined by independent methods). Before the code is applied, data for the cluster members (photometric brightness and colour) need to be prepared together with a file containing the list of clusters together with additional parameters. Examples of the data files are provided with the code. The code is applicable to Johnson (V,B-V), Gaia (G,BP-RP) and 2MASS (J,J-Ks) photometric systems.

The run-time of the code will depend on the total number of clusters in the included list, on the number of cluster members, and on the user input parameters. For example, if user inputs:

* delta_logAge = 0.2
* delta_Z = 0.005
* Nredd = 5
* Niter = 6

code will return results within 1-2 min for a typical open cluster. However, it can run for longer in the case of a larger cluster (for the included example of NGC 6791, the code returned results after 20 min). Furthermore, the specific run-time is also hardware-dependent (the code was tested on AMD Ryzen 3 PRO 4450U).

See Piecka & Paunzen (submitted) for a full description of the methods applied in the code.


## Requirements
In order to run the code, user must have installed Python 3 with numpy. The other libraries (matplotlib, time, os) are not required for the proper functionality of the code, but they provided additional information which may be useful.

The code was prepared for Windows 10, although it is also tested and usable on Linux. To make it work, only the lines depending on the os library need to be changed.
