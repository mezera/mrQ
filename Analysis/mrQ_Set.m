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
%   ~INPUTS~
%       mrQ:   The mrQ structure
%  saveflag:   The default is to save; to not save, enter 0
%     param:    {  SEE  }
%  varargin:    { BELOW }
%
% ~OUTPUTS~
%      mrQ:    The updated mrQ structure 
%
%    Parameter list and associated arguments:
%
%     ~PARAM~              ~VARARGIN~ 
%  1) alignflag
%  2) arrangerawflag
%  3) autoacpc
%  4) brakeaftert1            - [no longer in use]
%  5) brakeaftervisualization - [no longer in use]
%  6) channels
%  7) check                   - enter 1, if an interactive image is used to detect movement
%  8) clobber                 - enter "true" or "false" to redo the analysis
%  9) coilweights
% 10) complexflag
% 11) fieldstrength
% 12) freesurfer
% 13) inputdata_seir
% 14) inputdata_spgr
% 15) interp
% 16) lsq
% 17) lw
% 18) mmpervox
% 19) muti_coil_im
% 20) pdfit_method           - for M0-to-PD fit. 1= T1 single coil (default); 2= multi-coils; 3= both
% 21) permutation
% 22) polydeg                - the polynomial degree to fit to the Gain (default is 3)
% 23) proclus                - [no longer in use]
% 24) reset_raw_dir
% 25) rawdir                 - directory of raw image: a string 'path', as defined in mrQ_Create
% 26) refim                  - a path to an image to use as a refernce: a string 'path'
% 27) runfreesurfer
% 28) seir                   - list of SEIR scan numbers: a cell of strings, e.g. {'001' 002'}
% 29) siemens                - [no longer in use]
% 30) siemensdicom           - [no longer in use]
% 31) skip
% 32) spgr                   - list of SPGR scan numbers: a cell of strings, e.g. {'001' 002'}
% 33) sub                    - the name of subject: a string 'name'
% 34) sungrid                - enter 1, to use the SunGrid calls
% 35) useabs
%
% The following are checkers to determine whether or not a step has been
% performed. Inputting 1 means the step has been done, so don't redo it,
% and inputting 0 means the step hasn't been done, so do (or redo) it.
%     ~PARAM~              ~VARARGIN~ 
% 36) calm0_done             - M0 calculation         (1 or 0)
% 37) segmentation           - segmentation           (1 or 0)
% 38) seir_done              - SEIR T1 fit            (1 or 0)
% 39) spgr_coilweight_done   - SPGR coil weighting    (1 or 0)
% 40) spgr_init_done         - SPGR init step         (1 or 0)
% 41) spgr_pdbuild_done      - PD build               (1 or 0)
% 42) spgr_pdfit_done        - PD fit                 (1 or 0)
% 43) spgr_t1fit_done        - T1 fit SPGR            (1 or 0)
%
%
%
% Examples:
% mrQ = mrQ_Set(mrQ,,'raw','/biac2/wandell2/data/WMDevo/ambliopia/sub7/QuantativeImaging/20121102_3488')
% mrQ = mrQ_Set(mrQ,'sub','sub7');
% mrQ = mrQ_Set(mrQ,'SPGR',{ '0009' '0010' '0011' '0012' });
% mrQ = mrQ_Set(mrQ,'seir',{ '0004' '0005' '0006' '0007' });
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
% 2015
%
%


% Remove spaces and upper case
param = mrvParamFormat(param);
% This ensures that regardless of the user's input, mrQ_Set will always
% check the lower-case version of the string against the cases listed below

switch(param)
    
    case {'alignflag','align'}
        mrQ.alignFlag=varargin;
    case {'arrangerawflag'}
        mrQ.arrangeRawFlag=varargin;  
    case{ 'autoacpc'}
        mrQ.autoacpc = varargin;       
    case {'brakeaftert1' ,'t1brake'}
        mrQ.brakeAfterT1=varargin;         
    case {'brakeaftervisualization' ,'viewbrake'}
        mrQ.brakeAfterVisualization=varargin;
    case {'channels'}
        mrQ.channels=varargin;
    case {'check'}
        mrQ.check=varargin;        
    case {'clobber'}
        mrQ.clobber=varargin;        
    case {'coilweights', 'coilweight'}
        mrQ.coilWeights=varargin;        
    case {'complexflag','complex' }
        mrQ.complexFlag=varargin;        
    case{ 'field','stregth','fieldstrength'}
        mrQ.fieldstrength = varargin;
    case{ 'freesurfer'}
        mrQ.freesurfer=varargin;        
    case {'inputdata_seir'}
        mrQ.inputdata_seir=varargin;
    case {'inputdata_spgr'}
        mrQ.inputdata_spgr=varargin;        
    case {'interp'}
        mrQ.interp=varargin;        
    case {'lsq', 'lsqfit'}
        mrQ.lsq=varargin;        
    case{ 'lw','wl'}
        mrQ.LW= varargin; 
    case {'mmpervox', 'pixdim'}
        mrQ.mmPerVox=varargin;             
    case{ 'muti_coil_im'}
        mrQ.muti_coil_im=varargin;             
    case{'pdfit_method'}
      mrQ.PDfit_Method=varargin;               
    case {'permutation'}
       mrQ.permutation=varargin;            
    case {'polydeg' ,'poly','deg'}
       mrQ.PolyDeg=varargin;       
    case {'proclus' ,'proclass'}
        mrQ.proclus=varargin;            
    case{ 'reset_raw_dir'}
        mrQ.RawDir=fileparts(mrQ.name);              
    case {'rawdir','raw'}
        mrQ.RawDir=varargin;
    case {'refim','ref'}
        mrQ.refIm=varargin;
    case {'runfreesurfer' }
        mrQ.runfreesurfer=varargin;                    
    case {'seir_seriesnumbers','seir'}
        mrQ.SEIR_seriesNumbers=varargin;
    case {'siemens'}
        mrQ.siemens=varargin;
    case {'siemensdicom'}
        mrQ.SiemensDicom=varargin;    
    case {'skip'}
        mrQ.skip=varargin;        
    case {'spge_seriesNumbers','spgr'}
        mrQ.SPGR_seriesNumbers=varargin;
    case {'sub','name'}
        mrQ.sub=varargin;
    case {'sungrid' ,'grid'}
        mrQ.SunGrid=varargin;
    case {'useabs', 'magnitude' 'abs'}
        mrQ.useAbs=varargin;
     
        %Check: Has mrQ recorded it as done?
    case {'seir_done' }
        mrQ.SEIR_done=varargin;
    case {'spgr_init_done' }
        mrQ.SPGR_init_done=varargin;
    case {'spgr_coilweight_done'}
        mrQ.SPGR_coilWeight_done=varargin;
    case{ 'spgr_t1fit_done'}
        mrQ.SPGR_T1fit_done=varargin;
    case{ 'calm0_done'}
        mrQ.calM0_done=varargin;
    case{ 'spgr_pdfit_done'}
        mrQ.SPGR_PDfit_done=varargin;
    case{ 'spgr_pdbuild_done'}
        mrQ.SPGR_PDBuild_done=varargin;  
    case{ 'segmentation','seg'}
        mrQ.segmentation=varargin;   
      
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
