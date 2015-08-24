# mrQ manual (version 2) #
**mrQ** is a software package designed to calculate MR parameters (T1 and PD) using spoiled gradient echo scans (SPGR, FLASH). Using T1 and PD maps, mrQ performs the evaluation of macromolecular tissue volume (MTV) and the apparent volume of the interacting water protons (VIP) as well as the water-surface interaction rate (SIR). This is **version 2** of the mrQ software package.

The software and the tissue parameters are described in the following article:

>Mezer A, Yeatman JD, Stikov N, Kay KN, Cho NJ, Dougherty RF, Perry LM, Parvizi J, Hua LH, Butts-Pauly K and Wandell BA. Quantifying the local tissue volume and composition in individual brains with magnetic resonance imaging. Nature Medicine **19**, 1667-1672 (2013).
http://www.nature.com/nm/journal/v19/n12/full/nm.3390.html?WT.ec_id=NM-201312

and in:
> Mezer A, Rokem A, Hastie T and Wandell B. "Proton density mapping: Removing receive-inhomogeneity using multi-coil information and T1 regularization". ISMRM 23rd Annual Meeting & Exhibition. Toronto, Ontario, Canada. 2 June 2015.   

For the patent:

>Mezer, Aviv A., Robert F. Dougherty, and Brian A. Wandell. Methods for Detecting Abnormalities and Degenerative Processes in Soft Tissue Using Magnetic Resonance Imaging (MRI). The Board Of Trustees Of The Leland Stanford Junior University, assignee. Patent US9002428 B2. 7 Apr. 2015.

For an application of mrQ, see:

>Yeatman JD, Wandell BA and Mezer A. Lifespan maturation and degeneration of human brain white matter. Nature Communications **5** (2014).  http://www.nature.com/ncomms/2014/140917/ncomms5932/abs/ncomms5932.html

For more information, please contact: 

>Aviv Mezer: aviv.mezer(AT)elsc.huji.ac.il
>
>Shai Berman: shai.berman(AT)mail.huji.ac.il  
>
>Jonathan Bain: jonathan.bain(AT)mail.huji.ac.il



## Contents ##

- <a href=#versions>Versions</a>
    - <a href=#version-1>Version 1</a>
    - <a href=#version-2>Version 2</a> 
- <a href=#mrq-analysis-overview>mrQ analysis overview</a>
- <a href=#software-requirements>Software Requirements</a>
    - <a href=#required-third-party-software>Required third-party software</a>
    - <a href=#optional-software>Optional software</a>
    - <a href=#matlab-code>Matlab code</a>
- <a href=#mr-scanning->MR Scanning</a>
    - <a href=#spoiled-gradient-echo-scans-spgr-flash>Spoiled gradient echo scans (SPGR, FLASH)</a>
    - <a href=#epi-spin-echo-inversion-recovery-scans-b1-mapping>EPI spin-echo inversion recovery scans (B1 mapping)</a>
- <a href=#running-mrq>Running mrQ</a>
  - <a href=#overview>Overview</a>
  - <a href=#running-mrq-with-scitran-nifti-files>Running mrQ with SciTran NIfTI files</a>
  - <a href=#running-mrq-with-other-nifti-files>Running mrQ with other NIfTI files</a>
  - <a href=#example-directories>Example directories</a>
  - <a href=#visualization>Visualization</a>
  - <a href=#alignment>Alignment</a>
- <a href=#T1-fit-nonlinear-vs-weighted-linear-least-squares>T1 fit: nonlinear vs. weighted-linear least-squares</a>
- <a href=#parallel-computing>Parallel computing</a>
- <a href=#matlab-toolboxes>MATLAB Toolboxes</a>
- <a href=#forum>Forum</a>





### Versions ###
##### Version 1 #####
Version 1 (v.1) is the code that was used by Mezer et al. in their [Nature Medicine article](http://www.nature.com/nm/journal/v19/n12/full/nm.3390.html?WT.ec_id=NM-201312) (2013). The code is available at: https://github.com/mezera/mrQ/tree/v1.0

##### Version 2 #####
Version 2 (v.2) was released in autumn 2015 and was developed by the [Mezer lab](http://elsc.huji.ac.il/mezer/home) of the Hebrew University of Jerusalem, Israel, in collaboration with the [VISTA lab](https://vistalab.stanford.edu/) of Stanford University, USA. 

Some of the important changes are highlighted below: 
- A cleaner and more streamlined code, whose modular structure makes it easier to run and whose extensive commenting makes it easier to read.
- Only NIfTI files are accepted as input. DICOM files are no longer compatible with mrQ. We recommend the SciTran DICOM-to-NIfTI converter (see <a href=#optional-software>Optional Software</a>), which preserves key header information. If you use another converter, you will have to manually enter this information imto the mrQ structure (see <a href=#running-mrq-with-other-nifti-files>Running mrQ with other NIfTI files</a>).
- In past versions, setting up mrQ necessitated the use of three functions (mrQ_Create, mrQ_Set, mrQ_run). Now, these functions are all rolled into the brand-new mrQ_run. See <a href=#running-mrq>Running mrQ</a>.
- Though it will decrease the run time, the use of parallel computing is no longer required. The code is built to utilize SunGrid, though the default for SunGrid is “off”. See <a href=#parallel-computing>Parallel computing</a>.
- The default for the T1-M0 fit is now the weighted linear least squares method, replacing the nonlinear least squares method in past versions. This reduces the run time and is comparable to the nonlinear method in both precision and accuracy. See <a href=#T1-fit-non-linear-vs-weighted-least-squares>T1 fit: non-linear vs. weighted least-squares</a> for more information and citation.
- The B1 fit, which is performed using a one-parameter local regression and whose calculation is faster and more robust.
- Increased usage of parfor-loops, which decrease runtime through the use of parallel computation.
- For the M0-PD fit, a a local T1 regularization is used.
- Changes to the PD-CSF normalization yield more robust results.
- A nice synthetic T1-weighted image, which is good for segmentation.
- Multi-coil information can be advantageous, but is no longer required.

mrQ v.2 can be found at the link: _**[[[INSERT LINK HERE]]]**_. We recommend you use the most recent, up-to-date version of mrQ.

### mrQ analysis overview ###
- mrQ will use the mrQ structure you create and save it to the subject’s directory.
- New directories will be created, including directories for data and quantitative fits.
- Images will be registered to each other.
- SEIR-EPI T1 will be computed (low resolution) 
- SPGR T1, M0, B1 maps, and a synthetic T1-weighted image will be computed.
- T1-weighted and quantitative T1 images will be combined to segment the brain tissue. 
- PD and coil gain will be fitted from the M0 image.
- Biophysical models will be applied to calculate VIP and SIR maps.

### Software Requirements ###
##### Required third-party software #####
- MATLAB - http://www.mathworks.com/products/matlab/
  - Select MATLAB toolboxes are required for running mrQ. See <a href=#matlab-toolboxes>MATLAB Toolboxes</a> below.
- ANTs - http://stnava.github.io/ANTs/ 
- FSL - http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/
- SPM8 - http://www.fil.ion.ucl.ac.uk/spm/software/spm8/

##### Optional software #####
- FreeSurfer - http://surfer.nmr.mgh.harvard.edu/ 
    - *Earlier versions of the code relied on FreeSurfer, but it is no longer required. You can choose to use FreeSurfer at particular junctions in the code, e.g. Segmentation.* 
- Parallel computing environment (e.g. SunGrid Engine)
    - *Not required, but it will reduce running time.*
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
2. All scans should have a single TR.
3. Scans should have a minimal TE, about 2 msec. (Longer TE generates T2* weighting)
4. *Optional*: Save the multi-coil information. 
   - *GE*: Change the scanner default by editing the saveinter cv: saveinter=1.
   - *Siemens*: To save the individual coil information:  System &#10142; Miscellaneous &#10142; Save Uncombined.
   - *Philips*: If you know how to implement these settings, please email us or post in the forum. 
5. Scan with the same prescan parameters for all SPGR scans. 
   - *GE*: Scan the highest SNR image first (flip angle = 10 degrees). For the next scan, choose manual pre-scan and perform the scan without changing the pre-scan parameters.
   - *Siemens*: For the other flip angles: Right-click the sequence in the protocol &#10142; Properties and then select a Tab called Execution. There is an option called 'scan without further preparation' and you need to check that box.  That will cause that sequence to copy the prescan parameters from the previous one.
   - *Philips*: If you know how to implement these settings, please email us or post in the forum. 

##### EPI spin-echo inversion recovery scans (B1 mapping) #####

Low-resolution T1 maps are used to correct for the B1 bias. We will acquire data to fit unbiased T1 maps and correct the bias in the SPGR scans. The T1 fit is based on the <a href=#matlab-code>aforementioned</a> article by Barral et al. (2010).

1. Scan four SEIR-EPI readout scans with four different inversion times (50, 400, 1200 and 2400 msec).
2. Each scan needs to be acquired with slab inversion. 
   - *GE*: Change the scanner default by editing the a_gzrf0 cv: a_gzrf0=0
   - *Siemens*: For slab selective: Routine &#10142; Concatenations and set that parameter to the number of slices.  That way the new slice will not be acquired before the old one is finished.
   - *Philips*: If you know how to implement these settings, please email us or post in the forum. 
3. Use fat suppression. Fat suppression is recommended to be spatial-spectral to avoid any slice-selective imperfections.
   - *GE*: This is the default when slices are less than 4 mm thick.
   - *Siemens*: There is a similar setting called "water excitation", but check with Siemens for technical specifics.
   - *Philips*: If you know how to implement these settings, please email us or post in the forum. 

Alternatively, you can provide your own B1 map (NIfTI) if it is in SPGR space.

### Running mrQ ###
##### Overview #####
In v.1, running mrQ necessitated that the mrQ structure be created, set and executed. These steps were done using the mrQ functions mrQ_Create, mrQ_Set and mrQ_run, respectively. In v.2, these functions have all been incorporated into the brand-new **mrQ_run**, although mrQ_Create and mrQ_Set are still operational as standalone functions and can be called before executing mrQ_run.    

The mrQ_run code is modular, so it is easy to go step by step and see what is being computed when. 

If you have a NIfTI B1 map in SPGR space, you can input it into mrQ_run.

In past versions, changing a default parameter would be done in mrQ_Set, before running mrQ_run. Now, such changes can be done directly as mrQ_run input. After the first five inputs, enter a cell array of pairs (the parameter and its value), separated by commas. See examples below.

##### Running mrQ with SciTran NIfTI files #####
We recommend using SciTran to convert your DICOM files to NIfTI files (see <a href=#optional-software>Optional Software</a>), as it preserves the key header information for use in mrQ. 

If the SciTran-processed NIfTI files are located in the directory "dataDir" and the desired output folder is "outDir", you can run mrQ with the following code:

```matlab
mrQ_run(dataDir, outDir) 
```

If you have a B1 file (NIfTI) in SPGR space, you can input it into mrQ_run. If the B1 map is located at "B1file":
```matlab
mrQ_run(dataDir, outDir, [], [], B1file)
```

In the following syntax, the parameter "autoacpc" is being changed to 0 (default was 1):

```matlab
mrQ_run(dataDir, outDir, [],[], B1file, {'autoacpc', 0})
```

Enter as many parameters as you want to change, all inside one cell array. Write them in pairs, listing the parameter name followed by its value, with everything separated by commas:
```matlab
mrQ_run(dataDir, outDir, [], [], B1file, {'autoacpc', 0, 'sungrid', 1, 'wl', 0, 'refim', RefImageFile})

```

Alternatively, these can be performed in the mrQ_Set function, though they would have to be written as separate commands:
```matlab
mrQ = mrQ_Set('autoacpc', 0);
mrQ = mrQ_Set('sungrid', 1);
mrQ = mrQ_Set('wl', 0);
mrQ = mrQ_Set('refim', RefImageFile); 
```

##### Running mrQ with other NIfTI files #####
If you are using NIfTI files that were not generated by the SciTran code, you will need to manually enter key header information (using mrQ_Set) before you can run mrQ_run. 
- For the SPGR data, in addition to the location of the files and a unique name string for each, you will need each one's TR, TE, flip angle and field strength. 
- For the SEIR data, in addition to the location of the files and a unique name string for each, you will need each one's TR, TE and IT. 

If the NIfTI files are located in the directory "dataDir" and the desired output folder is "outDir", you can use the following code.

First, create a structure called "inputData_spgr". Set the required SPGR parameters: 

```matlab
%        A. Define the SPGR header info:
%
% dataDir is the location where the NIfTI are saved.
inputData_spgr.rawDir = dataDir;
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
```

Next, create a structure called "inputData_seir". Set the required SEIR parameters:
```matlab
%        B. Define the SEIR header info:
%
% dataDir is the location where the NIfTI are saved
inputData_seir.rawDir = dataDir;
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
```

Now we are ready to run mrQ:
```matlab
mrQ_run(dataDir, outDir, inputData_spgr, inputdata_seir,[])
```

If you have a B1 file (NIfTI) in SPGR space, you can input it into mrQ_run. If the B1 map is located at "B1file":
```matlab
mrQ_run(dataDir, outDir, inputData_spgr, inputdata_seir, B1file)
```

You can also change the parameters in this command:
```matlab
mrQ_run(dataDir, outDir, inputData_spgr, inputdata_seir, B1file, {'autoacpc', 0}) 
% or add as many parameters you want, as before, in one cell array and separated by commas
```

##### Example directories #####
The example directory and runscript for v.1 can be found at: http://purl.stanford.edu/qh816pc3429.

The example directory and runscript for v.2 can be found at: http://purl.stanford.edu/nn554zr6949.  
 - Note that the files in the v.2 directory are NIfTIs while those in the v.1 directory are DICOMs.

##### Visualization #####
To interactively watch the data after it's been aligned, define a "check" field and set it to 1. 
(It will activate the interaction in mrQ_initSPGR_ver2.m)
```matlab
mrQ_run(dataDir, outDir, [],[],[], {'check', 1}) 
```

##### Alignment #####
If a reference image is provided, the alignment in the SPGR section will be performed using that image. If no image is provided, mrQ will perform the AC-PC alignment automatically (default, 'acpc'=1), unless the user decides to do so manually (change to 'acpc'=0). Though the default is automatic, we recommend manual alignment when possible.

### T1 fit: nonlinear vs. weighted-linear least-squares ###
The most demanding computation in mrQ is the T1 fit. To avoid the long computing time of the nonlinear least squares method (which may take days on a single CPU for a whole brain with 1 mm<sup>3</sup> voxels), one can use the weighted-linear least-squares method as good alternative. The weighted linear method's accuracy and precision are comparable to those of the nonlinear method.

See: [Linear least-squares method for unbiased estimation of T1 from SPGR signals. Chang LC, Koay CG, Basser PJ, Pierpaoli C. Magn Reson Med. 2008 Aug;60(2):496-501. doi: 10.1002/mrm.21669.](http://www.ncbi.nlm.nih.gov/pubmed/18666108)

The default is to run the weighted-linear least-squares (and not the nonlinear least-squares). To reverse these settings (i.e. weighted-linear is off and nonlinear is on), use the following example code:
```matlab
mrQ_run(dataDir, outDir, [],[],[], {'wl', 0, 'lsq', 1})
```

### Parallel computing ###
Previous versions of the code relied on parallel computing for certain steps. However, v.2 of mrQ does not require parallel computing. If parallel computing such as the SunGrid engine is available, we recommend using it as a way of reducing the run time.

mrQ takes advantage of parallel computing during three steps within the analysis:

1. Calculating the transmit inhomogeneity for each voxel.
2. Calculating the T1 and M0 for each voxel.
3. Calculating the coil gain for different blocks in image space.

mrQ is written to take advantage of the SunGrid Engine for parallel computing, though its default is "off". Each user will need to change the specific calls to the grid, according to the parallel computing environment available. To turn SunGrid on, use the following example code:
```matlab
mrQ_run(dataDir, outDir, [],[],[], {'sungrid', 1}) 
```

### MATLAB Toolboxes ###
In addition to MATLAB, several MATLAB toolboxes are also required to run mrQ. The following nine toolboxes are required: (1) 'Bioinformatics Toolbox', (2) 'Image Processing Toolbox', (3) 'Optimization Toolbox', (4) 'Parallel Computing Toolbox', (5) 'Signal Processing Toolbox', (6) 'Simulink', (7) 'Statistical Parametric Mapping', (8) 'Statistics and Machine Learning Toolbox', (9) 'Symbolic Math Toolbox'.  

Certain toolboxes may only be required for optional segments in the code, though at the moment we cannot list exactly where each toolbox is used. 

Three notes on toolbox usage in mrQ v.2:
  -  In the function mrQ_boxScaleGlobLinear.m, the command "graphconncomp" is called. This is part of the Bioinformatics Toolbox.
  -  In the parallel computing (such as in mrQ_T1M0_LWFit.m), mrQ checks whether the Parallel Computing Toolbox is available. If so, it uses parfor-loops, as the "parfor" command is provided by this toolbox; if not, it uses for-loops. Though in these instances the Parallel Computing Toolbox is not required, it may be required at other junctions in the code.
  - In the function mrQ_T1M0_Fit.m, the function constructpolynomialmatrix3d.m is called, and this function uses the command "sym", which is provided by the Symbolic Math Toolbox. In this case, one can comment the usage of "sym" in constructpolynomialmatrix3d.m (lines 56-58) and proceed. However, note that the Symbolic Math Toolbox may be required at other junctions in the code.

Further updates about MATLAB toolbox usage will be announced on the forum or in future releases of mrQ.


### Forum ###
An interactive forum is available at _**[[[INSERT LINK HERE]]]**_. We hope that mrQ users will find it an informative site to discuss, share and troubleshoot for all topics relating to the mrQ software package.
