function mrQ = mrQ_Set(mrQ,param,varargin,saveflag)
% function mrQ = mrQ_Set(mrQ,param,varargin,saveflag)
%
% Setting data in the mrQ structure for the mrQ analysis.
%
% The mrQ structure stores all the mrQ computations.  mrQ_Set is the main
% function to add values, images, computations etc. to this structure.  As
% inputs, mrQ_Set accepts a mrQ structure ("mrQ"); a parameter ("param"), a
% string which is to be set; and a number of additional arguments and
% values that are specific to that parameter.  The arguments for each
% parameter are described below.  mrQ_Set always returns the modified mrQ
% structure and saves it, unless saveflag=0.  The structure is created
% using the default values from mrQ_Create or from the first call.
%
% See also:  mrQ_Create
%
%    Parameter list and associated arguments:
% saveflag         - the default is to save; to not save, enter 0
%
% ~PARAM~          ~VARARGIN~ (the key parameters are):
% rawdir               - directory of raw image: a string 'path', as defined in mrQ_Create
% seir                 - list of SEIR scan number: a cell of strings, e.g. {'001' 002'}
% SPGR                 - list of SPGR scan number: a cell of strings, e.g. {'001' 002'}
% sub                  - a name of subject: a string 'name'
% clobber              - enter "true" or "false" to redo the anlysis
% refim                - a path to an image to use as a refernce: a string 'path'
% check                - enter 1, if an interactive image is used to detect movement
% proclus              - enter 1, to use the proclus calls
% sungrid              - enter 1, to use the sungrid calls
% polydeg              - the polynomial degree to fit to the Gain (default is 3)
%
% brakeaftervisualization or viewbrake 
%                       - use that to start the data
% visual the images with the check flag and then stop the code from
% To reverse sewpesipic steps in the fir
%
% seir_done             - enter 1 to skip the SEIR T1 fit, or enter 0 to redo it
% spgr_init_done        - enter 1 to skip the SPGR init step, or enter 0 to redo it
% spgr_coilweight_done  - enter 1 to skip weighting, or enter 0 to redo it
% spgr_t1fit_done       - enter 1 to skip T1 fit SPGR, or enter 0 to redo it
% segmentation          - enter 1 to skip segmentation, or enter 0 to redo it
% calm0_done            - enter 1 to skip M0 calculation, or enter 0 to redo it
% spgr_pdfit_done       - enter 1 to skip PD fit, or enter 0 to redo it
% spgr_pdbuild_done     - enter 1 to skip PD build, or enter 0 to redo it
%
% ...                  - check code for others
%
% example:
% mrQ = mrQ_Set(mrQ,,'raw','/biac2/wandell2/data/WMDevo/ambliopia/sub7/QuantativeImaging/20121102_3488')
% mrQ = mrQ_Set(mrQ,'sub','sub7');
% mrQ = mrQ_Set(mrQ,'SPGR',{ '0009' '0010' '0011' '0012' });
% mrQ = mrQ_Set(mrQ,'seir',{ '0004' '0005' '0006' '0007' });
% mrQ = mrQ_Set(mrQ,'proclass',1);
%
%
% (C) Mezer lab, Hebrew University of Jerusalem, Israel
% 2015
%
%


% remove spaces and upper case
param = mrvParamFormat(param);

switch(param)
    
    
    case {'rawdir','raw'}
        mrQ.RawDir=varargin;
    case {'seir_seriesnumbers','seir'}
        mrQ.SEIR_seriesNumbers=varargin;
    case {'spge_seriesNumbers','spgr'}
        mrQ.SPGR_seriesNumbers=varargin;
    case {'sub','name'}
        mrQ.sub=varargin;
    case {'clobber'}
        mrQ.clobber=varargin;
    case {'arrangerawflag'}
        mrQ.arrangeRawFlag=varargin;
    case {'channels'}
        mrQ.channels=varargin;
    case {'complexflag','complex' }
        mrQ.complexFlag=varargin;
    case {'useabs', 'magnitude' 'abs'}
        mrQ.useAbs=varargin;
    case {'coilweights', 'coilweight'}
        mrQ.coilWeights=varargin;
    case {'alignflag','align'}
        mrQ.alignFlag=varargin;
    case {'interp'}
        mrQ.interp=varargin;
    case {'refim','ref'}
        mrQ.refIm=varargin;
    case {'mmpervox', 'pixdim'}
        mrQ.mmPerVox=varargin;
    case {'skip'}
        mrQ.skip=varargin;
    case {'check'}
        mrQ.check=varargin;
    case {'lsq', 'lsqfit'}
        mrQ.lsq=varargin;
    case {'runfreesurfer' }
        mrQ.runfreesurfer=varargin;
    case {'proclus' ,'proclass'}
        mrQ.proclus=varargin;
    case {'sungrid' ,'grid'}
        mrQ.SunGrid=varargin;
    case {'polydeg' ,'poly','deg'}
        mrQ.PolyDeg=varargin;
    case {'brakeaftervisualization' ,'viewbrake'}
        mrQ.brakeAfterVisualization=varargin;
    case {'brakeaftert1' ,'t1brake'}
        mrQ.brakeAfterT1=varargin;
    case {'siemens'}
        mrQ.siemens=varargin;
    case {'siemensdicom'}
        mrQ.SiemensDicom=varargin;
    case {'permutation'}
        mrQ.permutation=varargin;
    case {'inputdata_seir'}
        mrQ.inputdata_seir=varargin;
    case {'inputdata_spgr'}
        mrQ.inputdata_spgr=varargin;
    case {'seir_done' }
        mrQ.SEIR_done=varargin;
    case {'spgr_init_done' }
        mrQ.SPGR_init_done=varargin;
    case {'spgr_coilweight_done'}
        mrQ.SPGR_coilWeight_done=varargin;
    case{ 'spgr_t1fit_done'}
        mrQ.SPGR_T1fit_done=varargin;
    case{ 'segmentation','seg'}
        mrQ.segmentation=varargin;
    case{ 'calm0_done'}
        mrQ.calM0_done=varargin;
    case{ 'spgr_pdfit_done'}
        mrQ.SPGR_PDfit_done=varargin;
    case{ 'spgr_pdbuild_done'}
        mrQ.SPGR_PDBuild_done=varargin;
    case{ 'freesurfer'}
        mrQ.freesurfer=varargin;
    case{ 'reset_raw_dir'}
        mrQ.RawDir=fileparts(mrQ.name);
    case{ 'muti_coil_im'}
        mrQ.muti_coil_im=varargin;
    case{ 'autoacpc'}
        mrQ.autoacpc = varargin;
    case{ 'field','stregth','fieldstrength'}
        mrQ.fieldstrength = varargin;
    case{ 'wl'}
             mrQ.LW= varargin;  
    otherwise
        
        error('Unknown mrQ parameter!');
end


%save mrQ
if notDefined('saveflag')
    saveflag=1;
end

if saveflag==0
else
    save(mrQ.name,'mrQ');
end
