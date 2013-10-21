mrQ manual
===
mrQ is a software package designed to calculate MR parameters (T1 and PD) using spoiled gradient echo scans (SPGR, FLASH). mrQ allows the evaluation of macromolecule tissue fraction (MTV) and the volume of interaction proton (VIP) as well as the surface interaction ratio (SIR). 

The software and the tissue parameters are describe in the following article

Mezer A, Yeatman JD, Stikov N, Kay K, Cho NJ, Dougherty R, Perry LM, Parvizi J, Hua L, Butts-Pauly K, Wandell BA. Measuring within the voxel: brain macromolecular tissue volume in individual subjects. Nature Medicine, 2013 (in-press).

and the following patent application

Improved methods for detecting abnormalities in soft tissue using magnetic resonance imaging (MRI),USSN 61/437,587

For more information please contact

avivmezer@gmail.com



CONTENTS
====
- Software Requirements
- 3rd Party software
- Matlab code
- MR Scanning
- Spoiled gradient echo scans (SPGR,FLASH)
- EPI Spin echo inversion recovery scan (B1 mapping)
- Data organization
- Running mrQ
- Versions
- Parallel computing
- mrQ analysis overview
- Scanner dicoms types



Software Requirements
==
3rd Party software 
MATLAB  http://www.mathworks.com/products/matlab/ 
ANTS : http://stnava.github.io/ANTs/ 
FSL  http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/ 
and/or
FreeSrurfer http://surfer.nmr.mgh.harvard.edu/ 
Parallel computing environment (e.g., Sun Grid_Engine [Not required, but it will make a big difference for running time]).
Matlab code
mrQ makes use of other distributed matlab repositories:
mrQ - https://github.com/mezera/mrQ
Vistasoft  - https://github.com/vistalab/vistasoft
KNKUtils (from Kendrick Kay) - https://github.com/kendrickkay/knkutils
Joëlle Barral matlab code. A modified version of this code is part of the mrQ code. The original code can be found at: http://www-mrsrl.stanford.edu/~jbarral/t1map.html

MR Scanning 
==
Spoiled gradient echo scans (SPGR,FLASH)

1. 2-4 SPGR (not fast SPGR) scans with multiple flip angles recommended (e.g, 4, 10, 20, 30).  
2. All scans should have a single TR (note that a higher TR will increase the SNR).
3. Minimal TE (<=2ms)
4. Save the multi-coil information. To do this on GE scanners, change the scanner default by editing the saveinter cv: saveinter=1.
Scan with the same prescan parameters for all SPGR scans. To do this scan the highest SNR image first (flip angle = 10). For the next scan choose manual pre-scan and perform the scan without changing the pre-scan parameters.

EPI Spin echo inversion recovery scan (B1 mapping)

Low resolution T1 maps are used to correct for the B1 bias. We will scan data to fit unbiased T1 and correct the bias in the SPGR scans.

The T1 fit is based on Juelle Buarlle’s article: http://onlinelibrary.wiley.com/doi/10.1002/mrm.22497/abstract
http://www-mrsrl.stanford.edu/~jbarral/t1map.html 
Modified code is integrated within the mrQ software. 

Scan four SEIR - epi readout scans with four different inversion times (50, 400, 1200, 2400).
Each scan needs to be acquired with slab inversion on
GE scanner’s should change the scanner default by editing the a_gzrf0 cv: a_gzrf0=0
Use fat suppression. Fat suppression should be ??special-spectral?? to avoid any slice selective imperfections. This is the default with GE scanners when slices are less than 4mm thick.

Data organization
==
Follow these guidelines when organizing your data:

Data should be in a single directory - “DATA”.
Within DATA a dicoms directory is needed.
Within the dicoms directory each scan should be in a separate directory.
All SEIR dicoms need to be in the dicom directory.  
For each SPGR scan at list one dicom is needed so that the header can be read in by mrQ. 
All SPGR files need to be in nifti format within the DATA directory.

See http://purl.stanford.edu/qh816pc3429 for an example of directory organization.


Running mrQ 
==
To run mrQ a mrQ structure needs to be created, set and executed.
For example see ‘runScript’ in  http://purl.stanford.edu/qh816pc3429

Create a structure
mrQ=mrQ_Create(path)
Set mrQ field 
mrQ=mrQ_Set(mrQ,field name,field value)

For a given data set where SEIR scans are organized into 4 folders named  '0005' '0006' '0007' '0008' and SPGR scans are organized into 4 folders named '0009' '0010' '0011' '0012' the following can serve as an example script: 

% define the SEIR scans by the session’s 4 characters
mrQ=mrQ_Set(mrQ,'SEIR',{'0005' '0006' '0007' '0008'})

% define the SPGR scans by the session 4 characters
mrQ=mrQ_Set(mrQ,'SPGR',{'0009' '0010' '0011' '0012'})

% make a subject name
mrQ=mrQ_Set(mrQ,'sub','Examp')

% run:
mrQ_run(mrQ.name).

Versions
==
Version 1 (v.1) is the code to replicate that was used in Nature medicine mezer at. el. 2013 article: https://github.com/mezera/mrQ/tree/v1.0

We recommend you use the most recent, up to date of mrQ. The most active area of development is in the way the coil sensitivities are calculated. After V.1, later versions do not rely on Freesurfer any longer. An article describing those changes is in preparation. 

Parallel computing
==
mrQ takes advantage of parallel computing in three steps within analysis.
To calculate transmit inhomogeneity for each voxel.
To calculate T1 and M0 to each voxel
To calculate the coil gain for different bloc in image space.

mrQ is written to take advantage of the Sun grid parallel computing engine. Each user will need to change the specific calls to the grid according to the parallel computing environment available. One can turn off all those parallel jobs by editing the following setting when creating the mrQ structure:

mrQ=mrQ_Set(mrQ,sungrid’,0);
mrQ=mrQ_Set(mrQ,’proclus’,0);

If parallel computing not available to you please contact us. We are currently working on a general version of the code that does not rely on parallel computations. 

mrQ analysis overview
==
 mrQ will use the mrQ structure you create and save it to the subject’s directory.
 New directories will be created, including directories for data and quantitative fits.
 Images will be register to each other.
 SEIR-EPI T1 will be computed (low resultion) 
SPGR T1, M0, B1 maps, and a synthetic T1-weighted image, will be computed.
T1-weighted and quantitative T1 images will be combined to segment the brain tissue. 
PD and coil gain will be fit from the M0 image.
Biophysical model will be applay to calculate VIP and SIR maps.



Scanner dicoms types
==
The mrQ software was built around GE dicoms. It is possible that different vendors have different conventions in saving dicom information (e.g., header information, data ordering).

We are currently working on making the code compatible with different vendor’s dicom conventions. Please let us know if you experience any issues with reading dicoms when using the software. 

