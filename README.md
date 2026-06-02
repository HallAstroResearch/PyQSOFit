## PyQSOFit: A code to fit the spectrum of quasar  

### Getting Started

__See the [example](https://nbviewer.org/github/legolason/PyQSOFit/blob/master/example/example.ipynb) demo notebook for a quick start tutorial__

We provide a brief guide of the Python QSO fitting code (PyQSOFit) to measure spectral properties of SDSS quasars. The code was originally translated from Yue Shen's IDL code to Python. The package includes the main routine, Fe II templates, an input line-fitting parameter list, host galaxy templates, and dust reddening map to extract spectral measurements from the raw fits. Monte Carlo or MCMC estimation of the measurement uncertainties of the fitting results can be conducted with the same fitting code. 

The code takes an input spectrum (observed-frame wavelength, flux density and error arrays) and the redshift as input parameters, performs the fitting in the restframe, and outputs the best-fit parameters and quality-checking plots to the paths specified by the user. 

The code uses an input line-fitting parameter list to specify the fitting range and parameter constraints of the individual emission line components. An example of such a file is provided in the ``example.ipynb``. Within the code, the user can switch on/off components to fit to the pseudo-continuum. For example, for some objects the UV/optical Fe II emission cannot be well constrained and the user may want to exclude this component in the continuum fit. The code is highly flexible and can be modified to meet the specific needs of the user.

Use this code at your own risk, as we are not responsible for any errors. But if you report bugs, it will be greatly appreciated.


### Installation

Move to the folder in which you want the PyQSOFit source code folder to be created, then type:

git clone https://github.com/HallAstroResearch/PyQSOFit

You can edit the code in this folder, BUT SEE THE ALL-CAPS NOTE BELOW.  
Now install the program where python can find it (in your python path):

cd PyQSOFit

python -m pip install .

AFTER YOU CHANGE THE CODE YOU WILL NEED TO RERUN THE ABOVE TWO COMMANDS FOR THOSE CHANGES TO TAKE EFFECT.  

AND BE SURE TO MOVE THE EXAMPLE PYTHON NOTEBOOK TO A WORKING DIRECTORY ON YOUR OWN MACHINE BEFORE CHANGING AND RUNNING IT!  (SEE 'Your own python notebook copy' BELOW)  That will avoid it getting overwritten if you download a new version of the software from Github.


### Example Usage

Run the example0242.ipynb notebook in the j0242 folder:

jupyter notebook example0242.ipynb

The original example.ipynb notebook can still be found in the example folder, for reference.

The output files are:


Parameter file: <basename>_pp.dat 

Continuum fitting components described by 14 parameters (search for "Parameters used for continuum fitting" in PyQSOFit.py):

         pp[0]:     norm_factor for the MgII Fe_template

         pp[1]:     FWHM for the MgII Fe_template

         pp[2]:     small shift of wavelength for the MgII Fe template

         pp[3:5]:   same as pp[0:2] but for the Hbeta/Halpha Fe template

         pp[6]:     (PL_norm) norm_factor for continuum f_lambda = PL_norm * (lambda/3000.0)^{PL_slope}

         pp[7]:     (PL_slope) slope for the power-law continuum

         pp[8:10]:  norm, Te and Tau_e for the Balmer continuum below 3646 A

         pp[11:13]: polynomial for the continuum


Spectrum file (wave, flux, error): <basename>_PQF-RAW.dat 

These are the 'prereduced' values after removing bad pixels, dereddening, spectral trim, and smoothing. Probably should rename to "PQF-prered.dat".


Fit and Normalized Spectra file: <basename>_PQF_ASCII.dat 

The columns of PQF_ASCII.dat are:

1        2       3           4        5        6          7          8        9
rest_wav Mod_Con Mod_C+Lines fl/Mod_C er/Mod_C fl/Mod_C+L er/Mod_C+L fl-Mod_C fl-Mod_C+L

where Mod means Model, C means continuum and L means lines.

The last two columns are new, and are the flux minus the continuum fit and the flux minus the (continuum+lines) fit. 


### Your own python notebook copy

Copy the example0242.ipynb notebook so that you can edit it without being overwritten:

In the same folder in which the folder PyQSOFit is found, create a new folder. Let's say it's called MyFits (exact name doesn't matter).

Copy example0242.ipynb to that folder (MyFits, in this case) and rename the notebook (to avoid confusion).  Let's say it's now called myexample.ipynb (again, exact name doesn't matter).

In myexample.ipynb, change the sys.path.append lines in the first code block to :

sys.path.append('../PyQSOFit')
sys.path.append('../PyQSOFit/src')
sys.path.append('../PyQSOFit/src/pyqsofit')

Then, from the MyFits directory, run as usual:  jupyter notebook myexample.ipynb


### Hall Research Group Fork Changes

202505:
- Incorporated LMS (Lucas M. Seaton) changes to PyQSOFit.py.
- Edited example.ipynb to work with those changes.
- Added LMS_pyQSO_call.py to show how Lucas calls PyQSOFit.
- Added HostDecompV2.0.0.py and PyQSOFitV2.1.6.py purely for reference.
- Renumbered to 2.1.6.0; this fork's versions will have a 4th digit to distinguish them.

202506:
- Created example0242.ipynb

202507:
- Created option to use different Fe II template for J0242.
- Created version 2.1.6.1
- Enabled renaming of parameter file to something other than qsopar.fits
- Changed name of FeII optical and UV templates to be input parameters

To Do:
- investigate negative polynomials

### Cite

Preferred code citation: [Guo, H., Shen, Y., Wang, S. 2018, ascl:1809.008](https://ui.adsabs.harvard.edu/abs/2018ascl.soft09008G/abstract).

Please also cite: [Shen, Y. et al. 2019, ApJS, 241, 34S](https://ui.adsabs.harvard.edu/abs/2019ApJS..241...34S/abstract)

If using new host decompistion tools (`host_prior=True`), please cite: [Ren, W. et al. 2024](https://ui.adsabs.harvard.edu/abs/2024arXiv240617598R/abstract)
