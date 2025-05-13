#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 14:29:54 2025

@author: zulken
"""

"""The code below is how I call PyQSO in my overall code for J2318. It is not
'entirely' functional without the rest of the code for J2318."""

"---------------------------|=====================|---------------------------"
"   #####################   | PyQSOFit : Guo 2023 |   #####################   "
"---------------------------|=====================|---------------------------" 
if fitPQF==True:
    print('- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')   
    print('*** Fitting with PyQSOFit')
    "----------|===== S0 - Set up modules & saving convention =====|----------"  
    import os, sys, timeit, warnings
    sys.path.append('../')
    from pyqsofit.PyQSOFit import QSOFit
    from astropy.io import fits
    from astropy.table import Table
    warnings.filterwarnings("ignore")
    QSOFit.set_mpl_style()
    # Show the versions so we know what works
    import astropy, lmfit, pyqsofit, emcee
    # print(f'astropy version:    {astropy.__version__}')
    # print(f'lmfit version:      {lmfit.__version__}')
    # print(f'PyQSOFit Version:   {pyqsofit.__version__}')
    # print(f'emcee Version:      {emcee.__version__}')
    # print(f'PyQSOFit Directory: {pyqsofit.__path__}')
    ['/opt/anaconda3/lib/python3.11/site-packages/pyqsofit']
    FolderRoot = os.getcwd()+'/OUTPUT/'+objInfo['shortObjName']+'/'
    #Determination of current time:
    year   = time.localtime(time.time())[0]
    month  = time.localtime(time.time())[1]
    day    = time.localtime(time.time())[2]

    YYYY,MM,DD = str(year),str(month),str(day)
    if lmsT.MagOrd(month)==0: MM = '0'+str(month)
    if lmsT.MagOrd(day)==0: DD = '0'+str(day) 
    Date = YYYY+MM+DD+'/'

    Folder_PQF = FolderRoot+Date+'PyQSO_spectra/'
    if not os.path.exists(Folder_PQF): os.makedirs(Folder_PQF)
    
    "---------|===== Step 1: Set up the model input parameters =====|---------"    
    """ Firstly, run the script below to produce the line list file, qsopar.fits, 
    containing lines and their constraints, which will be needed in the following 
    fitting program. From this file, you can change some specific parameters to 
    suit your requirements, e.g., fitting range, line width, tie line center, tie 
    line sigma, etc. If you want to fit extra lines, please append it to 
    corresponding complex. Note that our line wavelength and sigma in the list are 
    in Ln scale, like Lnlambda, Lnsigma. """


    FilePath = '/Users/zulken/Desktop/PhD/Coding/InitialPlotting/PyQSOfit/PyQSOFit-v/example/data/'
    path_ex = Folder_PQF; #print(f'path_ex: {path_ex}')

    # create a header
    hdr0 = fits.Header()
    hdr0['Author'] = 'Lucas M. Seaton'
    primary_hdu = fits.PrimaryHDU(header=hdr0)
    ###############################################
    ######### Input 5
    ###############################################
    """ In this table, we specify the priors / initial conditions and boundaries 
    for the line fitting parameters. """
    
    #### NOTE: Only H\alpha should be shown for illustrative purposes
    
    line_priors = np.rec.array([
    #(𝜆_vac,name,𝜆min,𝜆max,linename,ngauss,𝐹_i,𝐹_min,𝐹_max,𝜎𝑖,𝜎min,𝜎max,voff,vindex,windex,findex,fvalue,vary)
     (6564.61, 'Ha', 6400, 6800, 'Ha_br', 2, 0.0, 0.0, 1e10, 5e-3, 0.004, 0.05, 0.015, 0, 0, 0, 0.05, 1),
     (6564.61, 'Ha', 6400, 6800, 'Ha_na', 1, 0.0, 0.0, 1e10, 1e-3, 5e-4, 0.00169, 0.01, 1, 1, 0, 0.002, 1),
     # (6549.85, 'Ha', 6400, 6800, 'NII6549', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 5e-3, 1, 1, 1, 0.001, 1),
     # (6585.28, 'Ha', 6400, 6800, 'NII6585', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 5e-3, 1, 1, 1, 0.003, 1),
     #(6718.29, 'Ha', 6400, 6800, 'SII6718', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 5e-3, 1, 1, 2, 0.001, 1),
     #(6732.67, 'Ha', 6400, 6800, 'SII6732', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 5e-3, 1, 1, 2, 0.001, 1),

     (4862.68, 'Hb', 4640, 5100, 'Hb_br', 2, 0.0, 0.0, 1e10, 5e-3, 0.004, 0.05, 0.01, 0, 0, 0, 0.01, 1),
     (4862.68, 'Hb', 4640, 5100, 'Hb_na', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 0.01, 1, 1, 0, 0.002, 1),
     (4960.30, 'Hb', 4640, 5100, 'OIII4959c', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 0.01, 1, 1, 0, 0.002, 1),
     (5008.24, 'Hb', 4640, 5100, 'OIII5007c', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 0.01, 1, 1, 0, 0.004, 1),
     (4960.30, 'Hb', 4640, 5100, 'OIII4959w',   1, 0.0, 0.0, 1e10, 3e-3, 2.3e-4, 0.004,  0.01,  2, 2, 0, 0.001, 1),
     (5008.24, 'Hb', 4640, 5100, 'OIII5007w',   1, 0.0, 0.0, 1e10, 3e-3, 2.3e-4, 0.004,  0.01,  2, 2, 0, 0.002, 1),
     #(4687.02, 'Hb', 4640, 5100, 'HeII4687_br', 1, 0.0, 0.0, 1e10, 5e-3, 0.004,  0.05,   0.005, 0, 0, 0, 0.001, 1),
     #(4687.02, 'Hb', 4640, 5100, 'HeII4687_na', 1, 0.0, 0.0, 1e10, 1e-3, 2.3e-4, 0.00169, 0.005, 1, 1, 0, 0.001, 1),

     #(3934.78, 'CaII', 3900, 3960, 'CaII3934' , 2, 0.0, 0.0, 1e10, 1e-3, 3.333e-4, 0.00169, 0.01, 99, 0, 0, -0.001, 1),

     #(3728.48, 'OII', 3650, 3800, 'OII3728', 1, 0.0, 0.0, 1e10, 1e-3, 3.333e-4, 0.00169, 0.01, 1, 1, 0, 0.001, 1),

     #(3426.84, 'NeV', 3380, 3480, 'NeV3426',    1, 0.0, 0.0, 1e10, 1e-3, 3.333e-4, 0.00169, 0.01, 0, 0, 0, 0.001, 1),
     #(3426.84, 'NeV', 3380, 3480, 'NeV3426_br', 1, 0.0, 0.0, 1e10, 5e-3, 0.0025,   0.02,   0.01, 0, 0, 0, 0.001, 1),

     (2798.75, 'MgII', 2700, 2900, 'MgII_br', 2, 0.0, 0.0, 1e10, 5e-3, 0.004, 0.05, 0.015, 0, 0, 0, 0.05, 1),
     (2798.75, 'MgII', 2700, 2900, 'MgII_na', 1, 0.0, 0.0, 1e10, 1e-3, 5e-4, 0.00169, 0.01, 1, 1, 0, 0.002, 1),

     #(1908.73, 'CIII', 1700, 1970, 'CIII_br', 2, 0.0, 0.0, 1e10, 5e-3, 0.004, 0.05, 0.015, 99, 0, 0, 0.01, 1),
     #(1908.73, 'CIII', 1700, 1970, 'CIII_na',   1, 0.0, 0.0, 1e10, 1e-3, 5e-4,  0.00169, 0.01,  1, 1, 0, 0.002, 1),
     #(1892.03, 'CIII', 1700, 1970, 'SiIII1892', 1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015,  0.003, 1, 1, 0, 0.005, 1),
     #(1857.40, 'CIII', 1700, 1970, 'AlIII1857', 1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015,  0.003, 1, 1, 0, 0.005, 1),
     #(1816.98, 'CIII', 1700, 1970, 'SiII1816',  1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015,  0.01,  1, 1, 0, 0.0002, 1),
     #(1786.7,  'CIII', 1700, 1970, 'FeII1787',  1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015,  0.01,  1, 1, 0, 0.0002, 1),
     #(1750.26, 'CIII', 1700, 1970, 'NIII1750',  1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015,  0.01,  1, 1, 0, 0.001, 1),
     #(1718.55, 'CIII', 1700, 1900, 'NIV1718',   1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015,  0.01,  1, 1, 0, 0.001, 1),

     #(1549.06, 'CIV', 1500, 1700, 'CIV_br', 2, 0.0, 0.0, 1e10, 5e-3, 0.004, 0.05, 0.015, 0, 0, 0, 0.05, 1),
     #(1549.06, 'CIV', 1500, 1700, 'CIV_na', 1, 0.0, 0.0, 1e10, 1e-3, 5e-4, 0.00169, 0.01, 1, 1, 0, 0.002, 1),
     #(1640.42, 'CIV', 1500, 1700, 'HeII1640',    1, 0.0, 0.0, 1e10, 1e-3, 5e-4,   0.00169, 0.008, 1, 1, 0, 0.002, 1),
     #(1663.48, 'CIV', 1500, 1700, 'OIII1663',    1, 0.0, 0.0, 1e10, 1e-3, 5e-4,   0.00169, 0.008, 1, 1, 0, 0.002, 1),
     #(1640.42, 'CIV', 1500, 1700, 'HeII1640_br', 1, 0.0, 0.0, 1e10, 5e-3, 0.0025, 0.02,   0.008, 1, 1, 0, 0.002, 1),
     #(1663.48, 'CIV', 1500, 1700, 'OIII1663_br', 1, 0.0, 0.0, 1e10, 5e-3, 0.0025, 0.02,   0.008, 1, 1, 0, 0.002, 1),

     #(1402.06, 'SiIV', 1290, 1450, 'SiIV_OIV1', 1, 0.0, 0.0, 1e10, 5e-3, 0.002, 0.05,  0.015, 1, 1, 0, 0.05, 1),
     #(1396.76, 'SiIV', 1290, 1450, 'SiIV_OIV2', 1, 0.0, 0.0, 1e10, 5e-3, 0.002, 0.05,  0.015, 1, 1, 0, 0.05, 1),
     #(1335.30, 'SiIV', 1290, 1450, 'CII1335',   1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015, 0.01,  1, 1, 0, 0.001, 1),
     #(1304.35, 'SiIV', 1290, 1450, 'OI1304',    1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.015, 0.01,  1, 1, 0, 0.001, 1),

     #(1215.67, 'Lya', 1150, 1290, 'Lya_br', 3, 0.0, 0.0, 1e10, 5e-3, 0.002, 0.05, 0.02, 0, 0, 0, 0.05, 1),
     #(1240.14, 'Lya', 1150, 1290, 'NV1240', 1, 0.0, 0.0, 1e10, 2e-3, 0.001, 0.01, 0.005, 0, 0, 0, 0.002, 1),
     # (1215.67, 'Lya', 1150, 1290, 'Lya_na', 1, 0.0, 0.0, 1e10, 1e-3, 5e-4, 0.00169, 0.01, 0, 0, 0, 0.002, 1),
     ],

        formats='float32,    a20,  float32, float32,      a20,  int32, float32, \
            float32, float32, float32, float32, float32, float32,   int32,  int32,  int32, float32, int32',
        names=' lambda, compname,   minwav,  maxwav, linename, ngauss,  inisca,  \
             minsca,  maxsca,  inisig,  minsig,  maxsig,    voff,  vindex, windex, findex,  fvalue,  vary')

    # Header
    hdr1 = fits.Header()
    hdr1['lambda'] = 'Vacuum Wavelength in Ang'
    hdr1['minwav'] = 'Lower complex fitting wavelength range'
    hdr1['maxwav'] = 'Upper complex fitting wavelength range'
    hdr1['ngauss'] = 'Number of Gaussians for the line'

    # Can be set to negative for absorption lines if you want
    hdr1['inisca'] = 'Initial guess of line scale [flux]'
    hdr1['minsca'] = 'Lower range of line scale [flux]'
    hdr1['maxsca'] = 'Upper range of line scale [flux]'

    hdr1['inisig'] = 'Initial guess of line sigma [lnlambda]'
    hdr1['minsig'] = 'Lower range of line sigma [lnlambda]'
    hdr1['maxsig'] = 'Upper range of line sigma [lnlambda]'

    hdr1['voff  '] = 'Limits on velocity offset from the central wavelength [lnlambda]'
    hdr1['vindex'] = 'Entries w/ same NONZERO vindex constrained to have same velocity'
    hdr1['windex'] = 'Entries w/ same NONZERO windex constrained to have same width'
    hdr1['findex'] = 'Entries w/ same NONZERO findex have constrained flux ratios'
    hdr1['fvalue'] = 'Relative scale factor for entries w/ same findex'

    hdr1['vary'] = 'Whether or not to vary the parameter (set to 0 to fix the line parameter to initial values)'

    # Save line info
    hdu1 = fits.BinTableHDU(data=line_priors, header=hdr1, name='line_priors')

    """ In this table, we specify the windows and priors / initial conditions and 
    boundaries for the continuum fitting parameters. """
    """Continuum Estimation Windows (RLF regions)"""
    conti_windows = np.rec.array(
        #[
        # (1150., 1170.), (1275., 1290.), (1350., 1360.), (1445., 1465.), # Continuum fitting 
        # (1690., 1705.), (1770., 1810.), (1970., 2400.), (2480., 2675.), # windows to avoid 
        # (2925., 3400.), (3775., 3832.), (4000., 4050.), (4200., 4230.), # emission line,
        # (4435., 4640.), (5100., 5535.), (6005., 6035.), (6110., 6250.), # etc.) [Å]
        # (6800., 7000.), (7160., 7180.), (7500., 7800.), (8050., 8150.), #
        
        # ],
        RLF,
        formats = 'float32,  float32',
        names =    'min,     max')

    hdu2 = fits.BinTableHDU(data=conti_windows, name='conti_windows')

    conti_priors = np.rec.array([
        ('Fe_uv_norm',  0.0,   0.0,   1e10,  1), # Normalization of the MgII Fe template [flux]
        ('Fe_uv_FWHM',  3000,  1200,  18000, 1), # FWHM of the MgII Fe template [Å]
        ('Fe_uv_shift', 0.0,   -0.01, 0.01,  1), # Wavelength shift of the MgII Fe template [ln(𝜆)]
        ('Fe_op_norm',  0.0,   0.0,   1e10,  1), # Normalization of the Hβ/H𝛼 Fe template [flux]
        ('Fe_op_FWHM',  3000,  1200,  18000, 1), # FWHM of the Hβ/H𝛼 Fe template [Å]
        ('Fe_op_shift', 0.0,   -0.01, 0.01,  1), # Wavelength shift of the Hβ/H𝛼 Fe template [ln(𝜆)]
        ('PL_norm',     1.0,   0.0,   1e10,  1), # Normalization of the power-law (PL) continuum f_𝜆=(𝜆/3000)^{-𝛼}
        ('PL_slope',    -1.5,  -5.0,  3.0,   1), # Slope of the power-law (PL) continuum
        ('Blamer_norm', 0.0,   0.0,   1e10,  1), # Normalization of the Balmer continuum at <3646 Å [flux]
        ('Balmer_Te',   15000, 10000, 50000, 1), # Te of the Balmer continuum at < 3646 Å [K?]
        ('Balmer_Tau',  0.5,   0.1,   2.0,   1), # Tau of the Balmer continuum at < 3646 Å
        ('conti_a_0',   0.0,   None,  None,  1), # 1st coefficient of the polynomial continuum
        ('conti_a_1',   0.0,   None,  None,  1), # 2nd coefficient of the polynomial continuum
        ('conti_a_2',   0.0,   None,  None,  1), # 3rd coefficient of the polynomial continuum
        # Note: The min/max bounds on the conti_a_0 coefficients are ignored by the code,
        # so they can be determined automatically for numerical stability.
        ],
        formats = 'a20,  float32, float32, float32, int32',
        names = 'parname, initial,   min,     max,     vary')

    hdr3 = fits.Header()
    hdr3['ini'] = 'Initial guess of line scale [flux]'
    hdr3['min'] = 'FWHM of the MgII Fe template'
    hdr3['max'] = 'Wavelength shift of the MgII Fe template'
    hdr3['vary'] = 'Whether or not to vary the parameter (set to 0 to fix the continuum parameter to initial values)'

    hdu3 = fits.BinTableHDU(data=conti_priors, header=hdr3, name='conti_priors')

    """
    In this table, we allow user to customize some key parameters in our result measurements.
    """
    measure_info = Table(
        [
            [[1350, 1450, 3000, 4200, 5100]],
            [[
                # [2240, 2650], 
                [4435, 4685],
            ]]
        ],
        names=([
            'cont_loc',
            'Fe_flux_range'
        ]),
        dtype=([
            'float32',
            'float32'
        ])
    )
    hdr4 = fits.Header()
    hdr4['cont_loc'] = 'The wavelength of continuum luminosity in results'
    hdr4['Fe_flux_range'] = 'Fe emission wavelength range calculated in results'

    hdu4 = fits.BinTableHDU(data=measure_info, header=hdr4, name='measure_info')
    hdu_list = fits.HDUList([primary_hdu, hdu1, hdu2, hdu3, hdu4])
    hdu_list.writeto(os.path.join(path_ex, 'qsopar.fits'), overwrite=True)

    # "Print the table:"
    # Table(line_priors)
    """NOTE: if you want to tie the line properties in different complex, you can 
    enlarge the complex range."""


    "-------------------|===== Step 2: Read the spectrum =====|-------------------"    
    """ Setup the paths and read in your spectrum. Our code is written under the 
    frame of SDSS spectral data format. Other data is also available as long as 
    they include wavelength, flux, error, and redshift, and make sure the 
    wavelength resolution is the same as SDSS spectrum (For SDSS the pixel scale is 
    1.e-4 in log space). """

    path_out = os.path.join(pyqsofit.__path__[0], '../', 'example/data/')
    # path_out = path_ex
    # print(f'path_out: {path_out}')

    # Requried
    key = 'UV2'
    data = cp.deepcopy(spectraObserved[key])
    lam = data[:,0]                         # OBS wavelength [Å]
    flux = data[:,1]                        # OBS flux [erg/s/cm^2/Å]
    err = data[:,2]                         # 1 sigma error
    z = objInfo['zem']                      # Redshift

    # Optional
    ra = objInfo['RA']                      # RA
    dec = objInfo['DEC']                    # DEC
    mjd = round(objInfo[key])               # SDSS MJD
    plateid = 104623                        # SDSS plate ID
    fiberid = 63050395075696130             # SDSS fiber ID
    RLF_PQF = cp.deepcopy(RLF)
    
    for key in spectra:
        if key =='UT1' or key == 'DT1':
            continue
        # key = 'UV2'                            #Force specific spectra
        mjd = round(objInfo[key])               # SDSS MJD
        data = cp.deepcopy(spectraObserved[key])
        lam = data[:,0]                         # OBS wavelength [Å]
        flux = data[:,1]                        # OBS flux [erg/s/cm^2/Å]
        err = data[:,2]                         # 1 sigma error

        "--------------------|===== Step 3: Fit the spectrum =====|-------------------"
        """ Use QSOFit to input the lam, flux, err, z, and other optimal parameters. 
        Use function Fit to perform the fitting. Default settings cannot meet all 
        needs. Please change settings for your own requirements. It depends on what 
        science you need. The following example set dereddening, host decomposition to 
        True.
    
        The broad_fwhm parameter can be adjusted depending on your definition (default 
        is 1200 km s−1). """
        # Prepare data
        q_mle = QSOFit(lam, flux, err, z, epoch=key, ra=ra, dec=dec, 
                       plateid=plateid, mjd=mjd, fiberid=fiberid, path=path_ex)
    
        # Double check the installation path with the PCA / Fe template files
        # print('install path:', q_mle.install_path)
    
        # Change it if you installed them somewhere else
        #q_mle.install_path = '...'
    
        start = timeit.default_timer()
        # Do the fitting
    
        q_mle.Fit(name=RMID,  # customize the name of given targets. Default: None -> plate-mjd-fiber
                  
                  # prepocessing parameters
                  nsmooth=1,            # do n-pixel smoothing to the raw input
                                        # flux and err spectra
                  and_mask=False,       # delete the and masked pixels
                  or_mask=False,        # delete the or masked pixels
                  reject_badpix=False,   # reject 10 most possible outliers by the 
                                        # test of pointDistGESD <- trims wavelength range
                  deredden=False,        # correct the Galactic extinction
                  wave_range=None,      # trim input wavelength
                  wave_mask=None,       # 2-D array, mask the given range(s)
    
                  # host decomposition parameters
                  #################################################################
                  # It is this parameter that is truncating the plot at ~3500AA
                  decompose_host=False,  # If True, the host galaxy-QSO 
                                        # decomposition will be applied (TRUE)
                  host_prior=True,      # If True, the code will adopt 
                                        # prior-informed  method to assist 
                                        # decomposition. Currently, only 'CZBIN1' 
                                        # and 'DZBIN1' model for QSO PCA are
                                        # available. And the model for galaxy must 
                                        # be PCA too.
                  host_prior_scale=0.2, # scale of prior panelty. Usually, 0.2 
                                        # works fine for SDSS spectra. Adjust it 
                                        # smaller if you find the prior affect the 
                                        # fitting results too much.
    
                  host_line_mask=True,  # If True, the line region of galaxy will be 
                                        # masked when subtracted from original spectra.
                  decomp_na_mask=True,  # If True, the narrow line region will be 
                                        # masked when perform decomposition
                  qso_type='CZBIN1',    # PCA template name for quasar
                  npca_qso=10,          # number of quasar templates
                  host_type='PCA',      # template name for galaxy
                  npca_gal=5,           # number of galaxy templates
                  
                  # continuum model fit parameters
                  Fe_uv_op=True,        # If True, fit continuum with UV and optical 
                                        # FeII template
                  poly=True,            # If True, fit continuum with the polynomial 
                                        # component to account for the dust reddening
                  #################################################################
                  # When BC=True the continuum is severely underestimated
                  BC=False,             # If True, fit continuum with Balmer continua 
                                        # from 1000 to 3646A (FALSE)
                  initial_guess=None,   # Initial parameters for continuum model, read 
                                        # the annotation of this function for detail
                  rej_abs_conti=False,  # If True, it will iterately reject 3 sigma 
                                        # outlier absorption pixels in the continuum
                  n_pix_min_conti=100,  # Minimum number of negative pixels for host 
                                        # continuuum fit to be rejected.
    
                  # emission line fit parameters
                  linefit=True,         # If True, the emission line will be fitted
                  rej_abs_line=False,   # If True, it will iterately reject 3 sigma 
                                        # outlier absorption pixels in the emission 
                                        # lines
    
                  # fitting method selection
                  MC=True,             # If True, do Monte Carlo resampling of the 
                                        # spectrum based on the input error array to 
                                        # produce the MC error array
                  MCMC=False,           # If True, do Markov Chain Monte Carlo sampling 
                                        # of the posterior probability densities to 
                                        # produce the error array. When True many WARNINGS are raised
                                        # WARNING:root:Too few points to create valid contours
                                        # The chain is shorter than 50 times the integrated autocorrelation time 
                                        # for 15 parameter(s). Use this estimate with caution and run a longer chain!
                  nsamp=200,            # The number of trials of the MC process (if 
                                        # MC=True) or number samples to run MCMC chain 
                                        # (if MCMC=True)
    
                  # advanced fitting parameters
                  param_file_name='qsopar.fits', # Name of the qso fitting parameter 
                                                 # FITS file.
                  nburn=20,                      # The number of burn-in samples to run 
                                                 # MCMC chain
                  nthin=10,                      # To set the MCMC chain returns every 
                                                 # n samples
                  epsilon_jitter=0.,             # Initial jitter for every initial 
                                                 # guess to avoid local minimum. (Under 
                                                 # test, not recommanded to change)
    
                  # customize the results
                  save_result=True,             # If True, all the fitting results will be 
                                                # saved to a fits file - Default: False
                  save_fits_name='FitResult',   # The output name of the result fits - Default: None
                  save_fits_path=Folder_PQF,    # The output path of the result fits - Default: path_out
                  plot_fig=True,                # If True, the fitting results will be 
                                                # plotted
                  save_fig=True,                # If True, the figure will be saved
                  plot_corner=True,             # Whether or not to plot the corner plot 
                                                # results if MCMC=True
    
                  # debugging mode
                  verbose=False,  # turn on (True) or off (False) debugging output
    
                  # sublevel parameters for figure plot and emcee
                  kwargs_plot={
                      # 'save_fig_path': '.',  # The output path of the figure
                      'save_fig_path': Folder_PQF
                      ,'broad_fwhm'   : 1200   # lower limit that code decide if a 
                                               # line component belongs to broad component [km/s]
                      ,'plot_residual': True  # Plotting of residuals
                      # ,'ylims' : [0.,150.]     # Limits of y-axis for ALL plots
                  },
                  kwargs_conti_emcee={},
                  kwargs_line_emcee={})
    
    
    """This end the PyQSO fitting. The following is data retrieval for various
    plotting and data dictionaries."""
    
    
    
    #Retrieve written ASCII. FUTURE: NEED to find way to ouput this as np.array,
    #                                rather than export it to .txt and reloading.
    ### wave_eval = np.linspace(np.min(self.wave) - 200, np.max(self.wave) + 200, 5000)
    WaveEval = False ### wave_eval - The wavelength range is different and thus so is the 
                     ###             flux and their errors.
    dp = Folder_PQF
    FL,PL = [],[]
    for path in os.listdir(dp):
        if path.endswith("PQF_RLF1Fix.dat"):
            FL.append(os.path.join(dp, path))
        if path.endswith("pp.txt"):
            PL.append(os.path.join(dp, path))
    conPyQSO = {}           # PyQSO continuum
    spectraPQF = {}         # normalized spectra
    PQFvalues = {}          # Best Fit Parameters
    for file in FL:
        ep = file[-19:-16]; #print(ep)
        
        data = cp.deepcopy(spectra[ep])
        lam,flux,err = data[:,0],data[:,1],data[:,2]
        
        conPyQSO[ep] = np.genfromtxt(file, usecols=(0,1,2)) 
        #NOTE: There is a problem where len(conPyQSO[ep]) is slightly smaller than len(flux)
        if len(conPyQSO[ep]) != len(flux): 
            pyqsofit = conPyQSO[ep]

        else:
            pyqsofit = np.zeros(np.shape(data))
            pyqsofit[:,0] = lam
            pyqsofit[:,1] = flux/conPyQSO[ep][:,1]
            pyqsofit[:,2] = err/conPyQSO[ep][:,1]
            for p,px in enumerate(lam): #Error propagation for PQF normalized flux uncertainties
                pyqsofit[p,2] = pyqsofit[p,1] * np.sqrt((err[p]/flux[p])**2 + (conPyQSO[ep][p,2]/conPyQSO[ep][p,1])**2)
        spectraPQF[ep] = pyqsofit
    for file in PL:
        ep = file[-10:-7]
        PQFvalues[ep] = np.genfromtxt(file, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13))
    ###########################################################################
    #### Performing Goodness of Fit Test
    ###########################################################################
    GFit_PQF,NPQF,CPQF,rCPQF,NPQF = {},[],[],[],[]
    CPQFBW =[['',0],['',np.inf],['',0],['',np.inf]] #Chi2 power law Best/Worst
    M_PQF = len(PQFvalues[ep])
    for ep in spectra:
        if ep == 'DT1' or ep == 'UT1':
            continue
        if WaveEval == True:
            continue
        GFit_PQF[ep] = lmsT.GoodFit2023LMS(spectra[ep], conPyQSO[ep][:,1], M_PQF, RLF_PQF, objInfo)
        CPQF.append(GFit_PQF[ep][0]); rCPQF.append(GFit_PQF[ep][1]); NPQF.append(GFit_PQF[ep][2])
        if GFit_PQF[ep][0]>CPQFBW[0][1]: CPQFBW[0] = [ep, GFit_PQF[ep][0]] #MAX chi2
        if GFit_PQF[ep][0]<CPQFBW[1][1]: CPQFBW[1] = [ep, GFit_PQF[ep][0]] #MIN chi2
        if GFit_PQF[ep][1]>CPQFBW[2][1]: CPQFBW[2] = [ep, GFit_PQF[ep][1]] #MAX reduced chi2   
        if GFit_PQF[ep][1]<CPQFBW[3][1]: CPQFBW[3] = [ep, GFit_PQF[ep][1]] #MIN reduced chi2
    nuPQF = round(np.mean(NPQF)) - M_PQF
    ###########################################################################
    #### Recording PyQSO data into FITinfo dictionary
    ###########################################################################
    # LPQFEQ = 'Continuum: PL + F_poly_conti'
    # LPQFEQ = r'$F(\lambda)=(\lambda/3000)^{-B}+FeUV+FeOP+\beta+poly$'
    LPQFEQ = r'$F(\lambda)=A(\lambda/3000)^{-B}+U+O+\beta+C$'
    FITinfo['model'].append([['PyQSOFit','PQF'], LPQFEQ, 'green',
                               [r'$U_n$', r'$U_{2\gamma}$', r'$U_\lambda$', r'$O_n$',
                                r'$O_{2\gamma}$', r'$O_\lambda$', r'$A$', r'$B$',
                                r'$\beta_n$', r'$\beta_T$', r'$\beta_\tau$',
                                r'$C_0$', r'$C_1$', r'$C_2$'],
                              [Chi,chi,r'$N$']])
    if WaveEval == True:
        CPQFBW = CBRKBW         #Adopt RBP2 values for chi2
        GFit_PQF = GFit_BRK     #Adopt RBP2 values for chi2
    FITinfo['MaxChi'].append(CPQFBW[0]); FITinfo['MinChi'].append(CPQFBW[1])
    for ep in spectra:
        if ep == 'DT1' or ep == 'UT1':
            continue
        if ep in FITinfo: 
            FITinfo[ep].append([PQFvalues[ep], GFit_PQF[ep], conPyQSO[ep][:,1], spectraPQF[ep]])
        else: 
            FITinfo[ep] = [[PQFvalues[ep], GFit_PQF[ep], conPyQSO[ep][:,1], spectraPQF[ep]]]
    ###########################################################################         
    ### Recording into NEW nested dictionaries, FITinfo2
    ###########################################################################
    model = 'PQF'
    FITinfo2['model'][model] = {}
    FITinfo2['model'][model]['name'] = 'PyQSO'
    FITinfo2['model'][model]['eqn'] = LPQFEQ
    FITinfo2['model'][model]['colour'] = 'green'
    FITinfo2['model'][model]['fit_param'] = [r'$U_n$', r'$U_{2\gamma}$', r'$U_\lambda$', r'$O_n$',
                                              r'$O_{2\gamma}$', r'$O_\lambda$', r'$A$', r'$B$',
                                              r'$\beta_n$', r'$\beta_T$', r'$\beta_\tau$',
                                              r'$C_0$', r'$C_1$', r'$C_2$']
    FITinfo2['model'][model]['good_fit'] = [Chi, chi, r'$N$']
    FITinfo2['model'][model]['RLF'] = RLF_PQF
    FITinfo2['model'][model]['FeII_RLF'] = FeII_RLF
    ### Record Fitting parameters, Chi2, continuum, and normalization for each epoch
    for ep in spectra:
        FITinfo2[ep][model] = {}
        FITinfo2[ep][model]['fit_param'] = PQFvalues[ep]
        FITinfo2[ep][model]['good_fit'] = GFit_PQF[ep]
        FITinfo2[ep][model]['continuum'] = conPyQSO[ep]
        FITinfo2[ep][model]['normalized'] = spectraPQF[ep]
    ### Recording Maxima and Minima Chi2 values
    record = 0
    for f,fit in enumerate(FITinfo2['model']):
        if fit == model and 0<=f<len(FITinfo2['MaxChi']):
            FITinfo2['MaxChi'][f] = CPQFBW[0]
            FITinfo2['MinChi'][f] = CPQFBW[1]
            record += 1
            break
    if record == 0:
        FITinfo2['MaxChi'].append(CPQFBW[0])
        FITinfo2['MinChi'].append(CPQFBW[1])
    
    end = timeit.default_timer()
    print(f'Fitting finished in {np.round(end - start, 1)}s')
   
    """
    The gray shaded bars at the top are the continuum windows used in the fitting.
    Now you are already done with the QSO fitting part! """