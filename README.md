# mrQ manual (version 2) #
mrQ is a software package designed to calculate MR parameters (T1 and PD) using spoiled gradient echo scans (SPGR, FLASH). mrQ enables the evaluation of macromolecular tissue volume (MTV) and the apparent volume of the interacting water protons (VIP) as well as the water-surface interaction rate (SIR). 

The software and the tissue parameters are described in the following article:

>Mezer A, Yeatman JD, Stikov N, Kay KN, Cho NJ, Dougherty RF, Perry LM, Parvizi J, Hua LH, Butts-Pauly K and Wandell BA. Quantifying the local tissue volume and composition in individual brains with magnetic resonance imaging. Nature Medicine **19**, 1667-1672 (2013).
http://www.nature.com/nm/journal/v19/n12/full/nm.3390.html?WT.ec_id=NM-201312

and the following patent application:

>Improved methods for detecting abnormalities in soft tissue using magnetic resonance imaging (MRI), USSN 61/437,587

For an application of mrQ, see:

>Yeatman JD, Wandell BA and Mezer A. Lifespan maturation and degeneration of human brain white matter. Nature Communications **5** (2014).  http://www.nature.com/ncomms/2014/140917/ncomms5932/abs/ncomms5932.html

For more information, please contact: 

>Aviv Mezer, avivmezer@gmail.com



## Contents ##

- <a href=#versions>Versions</a>
    - <a href=#version-1>Version 1</a>
    - <a href=#version-2>Version 2</a> 
- <a href=#mrq-analysis-overview>mrQ analysis overview</a>
- <a href=#software-requirements>Software Requirements</a>
    - <a href=#required-third-party-software>Required third-party software</a>
    - <a href=#optional-software>Optional Software</a>
    - <a href=#matlab-code>Matlab code</a>
- <a href=#mr-scanning->MR Scanning</a>
    - <a href=#spoiled-gradient-echo-scans-spgrflash>Spoiled gradient echo scans (SPGR, FLASH)</a>
    - <a href=#epi-spin-echo-inversion-recovery-scans-b1-mapping>EPI spin-echo inversion recovery scans (B1 mapping)</a>
- <a href=#data-organization>Data Organization</a>
- <a href=#running-mrq>Running mrQ</a>
  - <a href=#overview>Overview</a>
  - <a href=#running-mrq-with-scitran-nifti-files>Running mrQ with SciTran NIfTI files</a>
  - <a href=#running-mrq-with-other-nifti-files>Running mrQ with other NIfTI files</a>
- <a href=#T1-fit-non-linear-vs-weighted-least-squares>T1 fit: non-linear vs. weighted least-squares</a>
- <a href=#parallel-computing>Parallel computing</a>


Running mrQ with SciTran NIfTI files


### Versions ###
##### Version 1 #####
Version 1 (v.1) is the code that was used by Mezer et al. in their Nature Medicine (2013) article. It is available at: https://github.com/mezera/mrQ/tree/v1.0

##### Version 2 #####
Version 2 (v.2) was released in autumn 2015 and was developed by the Mezer lab of the Hebrew University of Jerusalem, Israel. Some of the important changes are highlighted below: 
- A cleaner and more streamlined code, whose modular structure makes it easier to run and whose extensive commenting makes it easier to read.
- Only NIfTI files are accepted as input. DICOM files are no longer compatible with mrQ. We recommend the SciTran DICOM-to-NIfTI converter (see <a href=#optional-software>Optional Software</a>), which preserves key header information. If you use another converter, you will have to manually enter this information imto the mrQ structure (see <a href=#running-mrq-with-other-nifti-files>Running mrQ with other NIfTI files</a>).
- The code runs without the need of parallel computing (though it will go much faster with parallel computing). The code is built to utilize SunGrid, though the default for SunGrid is “off”. See <a href=#parallel-computing>Parallel computing</a>.
- More new defaults include using the weighted linear least squares for the T1-M0 fit; using local regression for the B1 fit; and using a local T1 regularization for the M0-PD fit.
- The B1 fit, which is performed using a one-parameter local regression.
- Changes to the PD-CSF normalization, yielding more robust results.
- A nice synthetic T1-weighted image, which is good for segmentation.

mrQ v.2 can be found at the link: _**[[[INSERT LINK HERE]]]]**_. We recommend you use the most recent, up-to-date version of mrQ.

### mrQ analysis overview ###
- mrQ will use the mrQ structure you create and save it to the subject’s directory.
- New directories will be created, including directories for data and quantitative fits.
- Images will be registered to each other.
- SEIR-EPI T1 will be computed (low resultion) 
- SPGR T1, M0, B1 maps, and a synthetic T1-weighted image will be computed.
- T1-weighted and quantitative T1 images will be combined to segment the brain tissue. 
- PD and coil gain will be fitted from the M0 image.
- Biophysical models will be applied to calculate VIP and SIR maps.

### Software Requirements ###
##### Required third-party software #####
- MATLAB - http://www.mathworks.com/products/matlab/ 
- ANTs - http://stnava.github.io/ANTs/ 
- FSL - http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/
- SPM8 - http://www.fil.ion.ucl.ac.uk/spm/software/spm8/

##### Optional software #####
- FreeSurfer - http://surfer.nmr.mgh.harvard.edu/ 
    - *Earlier versions of the code relied on FreeSurfer, but it is no longer required. You can choose to use FreeSurfer at particular junctions in the code, e.g. Segmentation.* 
- Parallel computing environment (e.g. SunGrid Engine)
    - *Not required, but it will make a big difference for running time.*
- A DICOM-to-NIfTI converter
    - *We recommend using the code that was developed by [Scientific Transparency](http://scitran.github.io/) (SciTran), a project at Stanford University. This code preserves key header information that is required by mrQ. To get the code, you will need to use [Docker](http://docs.docker.com/windows/started/). The code is short and is available at: https://hub.docker.com/r/vistalab/nimsdata/*
    - *You can use other types of NIfTI files, but you will need to manually enter the missing header information.  See* <a href=#running-mrq-with-other-nifti-files>Running mrQ with other NIfTI files</a>. 


##### Matlab code #####
mrQ requires the following openly distributed code repositories:

- mrQ - https://github.com/mezera/mrQ
    - The mrQ repository contains a modified version of Joëlle Barral's matlab code. The original code can be found at http://www-mrsrl.stanford.edu/~jbarral/t1map.html and was originally published in 
[Barral JK, Gudmundson E, Stikov N, Etezadi-Amoli M, Stoica P, and Nishimura DG (2010). A robust methodology for in vivo T1 mapping. Magnetic Resonance in Medicine, Oct. 60(4), 1057-1067. doi:10.1002/mrm.22497](http://onlinelibrary.wiley.com/doi/10.1002/mrm.22497/abstract)
- Vistasoft  - https://github.com/vistalab/vistasoft
- KNKUtils (from Kendrick Kay) - https://github.com/kendrickkay/knkutils


### MR Scanning ###
##### Spoiled gradient echo scans (SPGR, FLASH) #####

1. 2-4 SPGR (not fast SPGR) scans with multiple flip angles recommended (e.g., 4, 10, 20 and 30 degrees).  
2. All scans should have a single TR (note that a higher TR will increase the SNR).
3. Minimal TE (less than or equal to 2 msec)
4. Save the multi-coil information. To do this on GE scanners, change the scanner default by editing the saveinter cv: saveinter=1.
5. Scan with the same prescan parameters for all SPGR scans. To do this, scan the highest SNR image first (flip angle = 10). For the next scan, choose manual pre-scan and perform the scan without changing the pre-scan parameters.

##### EPI spin-echo inversion recovery scans (B1 mapping) #####

Low-resolution T1 maps are used to correct for the B1 bias. We will acquire data to fit unbiased T1 maps and correct the bias in the SPGR scans. The T1 fit is based on the <a href=#matlab-code>aforementioned</a> article by Barral et al. (2010).

1. Scan four SEIR-EPI readout scans with four different inversion times (50, 400, 1200 and 2400 msec).
2. Each scan needs to be acquired with slab inversion. For GE scanners, change the scanner default by editing the a_gzrf0 cv: a_gzrf0=0
3. Use fat suppression. Fat suppression should be spatial-spectral to avoid any slice-selective imperfections. Note: This is the default for GE scanners when slices are less than 4 mm thick.

### Data Organization ###
#### Follow these guidelines when organizing your data: ####

- Data should be in a single directory, “DATA”.
- Within the DATA directory, a DICOMs directory is needed.
- Within the DICOMs directory, each scan should be in a separate directory.
- All SEIR DICOMs need to be in the DICOM directory.  
- A single DICOM is needed for each scan in that scan's directory so that the header can be read by mrQ. 
- All SPGR files need to be in NIfTI format within the DATA directory.

See http://purl.stanford.edu/qh816pc3429 for an example of directory organization.


### Running mrQ ###
##### Overview #####
To run mrQ, a mrQ structure (a .mat file) needs to be created, set and executed. These steps are done using the mrQ functions **mrQ_Create**, **mrQ_Set** and **mrQ_run_Ver2**, respectively. 
- The function **mrQ_Create** takes the data files' location ("dataDir"), creates an output directory ("outDir") and initializes the mrQ structure.
  ```matlab           
  mrQ = mrQ_Create(dataDir, [], outDir)
  ```       
- The function **mrQ_Set** allows the user to set a number of parameters, such as the desired method of the PD fit or whether the SunGrid is available. If your NIfTI files were not created using the SciTran code, you will need to enter the missing information here.
    ```matlab
    mrQ = mrQ_Set(mrQ,parameter_name,paramter_value)
    % for example:
    mrQ = mrQ_Set(mrQ, pdfit_method, 2)
     ```
- The function **mrQ_run_Ver2** will run mrQ. If the NIfTI files were created using the SciTran code, the location of the NIfTI files is sufficient to run the code with the default parameters. 
   ```matlab
   mrQ_run(mrQ.name)
   ```
   The mrQ_run_Ver2 code is modular, so it is easy to go step by step and see what is being computed when. 
   
   An added functionality of the v.2 run code is parameters can be set directly as input of mrQ_run_Ver2, and not only through mrQ_Set. For example:
   ```matlab
   mrQ_run(mrQ.name, pdfit_method, 2)
   ```

For an example of this structure, see ‘runScript’ at: http://purl.stanford.edu/qh816pc3429

#### Running mrQ with SciTran NIfTI files ####
We recommend using SciTran to convert your DICOM files to NIfTI files (see <a href=#optional-software>Optional Software</a>), as it preserves the key header information for use in mrQ. Once the NIfTI files are in their directory, you can run mrQ. For example:

```matlab
% Create mrQ structure and define the datadir where the NIfTI are saved, 
% and outdir where mrQ output will be saved.
mrQ = mrQ_Create(dataDir,[],outDir);
%
% One can set many different fit properties via mrQ_set.m
% 
% Get the image and hdr info from the NIfTI:
mrQ = mrQ_arrangeData_nimsfs(mrQ);
%
% Run it:
mrQ_run(mrQ.name)
```

#### Running mrQ with other NIfTI files ####
If you are using NIfTI files that were not generated by the SciTran code, you will need to manually enter key header information using mrQ_Set. 
- For the SPGR data, in addition to the location of the files and a unique name string for each, you will need each one's TR, TE, flip angle and field strength. 
- For the SEIR data, in addition to the location of the files and a unique name string for each, you will need each one's TR, TE and IT. 

Example:
```matlab
% Create mrQ structure and define the datadir where the NIfTI are saved,
% and outdir where mrQ output will be saved.
mrQ = mrQ_Create(dataDir,[],outDir);
%
% Set different properties via mrQ_Set, if desired.
%
%        A. Define the SPGR header info:
%
% mrQ.RawDir is the location where the NIfTI are saved.
inputData_spgr.rawDir = mrQ.RawDir;
%
% A list of NIfTI names. (A unique string from the names is enough)
inputData_spgr.name = {'0009' '0010' '0011' '0012'};
%
% The TR of each NIfTI in the list (in msec)
inputData_spgr.TR = [12 12 12 12];
%
% The TE of each NIfTI in the list (in msec)
inputData_spgr.TE = [2.27 2.27 2.27 2.27];
%
% The flip angle of each NIfTI in the list (in degrees)
inputData_spgr.flipAngle = [4 10 20 30];
%
% The field strength of each NIfTI in the list (in Teslas)
inputData_spgr.fieldStrength = [3 3 3 3];
%
%        B. Define the SEIR header info:
%
% mrQ.RawDir is the location where the NIfTI are saved
inputData_seir.rawDir = mrQ.RawDir;
%
% A list of NIfTI names.  (A unique string from the names is enough)
inputData_seir.name = {'0005'  '0006'  '0007'  '0008'};
%
% The TR of each NIfTI in the list (in msec)
inputData_seir.TR = [3000 3000 3000 3000];
%
% The TE of each NIfTI in the list (in msec)
inputData_seir.TE = [49 49 49 49];
%
% The inversion time of each NIfTI in the list (in msec)
inputData_seir.IT = [50 400 1200 2400];
%
% Add the NIfTI info to the mrQ structure:
mrQ = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir);
%
% Run it:
mrQ_run(mrQ.name)
```

##### Visualization #####
To interactively watch the data after it's been aligned, define a *check* field and set it to 1. 
(It will activate the interaction in mrQ_initSPGR_ver2.m)

##### Alignment #####
The default alignment of the images is an automatic AC-PC alignment. 
It can be semi-manual or non-AC-PC. Settings can be changed when calling mrQ_run or by changing the settings later with mrQ_Set.

### T1 fit: non-linear vs. weighted least-squares ###
The most demanding computation in mrQ is the T1 fit. To avoid the long computing time of the non-linear method (which may take days on a single CPU for a whole brain with 1 mm<sup>3</sup> voxels), one can use the weighted-least square method as good alternative.

See: [Linear least-squares method for unbiased estimation of T1 from SPGR signals. Chang LC, Koay CG, Basser PJ, Pierpaoli C. Magn Reson Med. 2008 Aug;60(2):496-501. doi: 10.1002/mrm.21669.](http://www.ncbi.nlm.nih.gov/pubmed/18666108)

To run the weighted least-squares:
```matlab
mrQ = mrQ_Set(mrQ,'wl',1);
mrQ = mrQ_Set(mrQ,'lsq',0);
```

### Parallel computing ###
Previous versions of the code relied on parallel computing for certain steps. However, mrQ does not require parallel computing in v.2 and later. If parallel computing such as the SunGrid engine is available, we recommend using it as a way of greatly reducing the run time.

mrQ takes advantage of parallel computing during three steps within the analysis:
1. Calculating the transmit inhomogeneity for each voxel.
2. Calculating the T1 and M0 for each voxel.
3. Calculating the coil gain for different blocks in image space.

mrQ is written to take advantage of the SunGrid Engine for parallel computing, though its default is "off". Each user will need to change the specific calls to the grid, according to the parallel computing environment available. The user can turn off or on all of those parallel jobs when creating the mrQ structure:
```matlab
mrQ=mrQ_Set(mrQ,'sungrid',0); % off, default
% or
mrQ=mrQ_Set(mrQ,'sungrid',1); % on
```

