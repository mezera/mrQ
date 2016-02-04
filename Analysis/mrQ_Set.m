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
%  1) alignflag               - Enter to 1 if you want to align the SEIR 
%                                   slices, 0 (default) if not. If the
%                                   SEIR data has more than a few slices,
%                                   it is a good idea to try to align them.
%  2) arrangerawflag          - [no longer in use]
%  3) autoacpc                - Enter 1 (default) if you want mrQ to 
%                                   perform the AC-PC alignment 
%                                   automatically, or enter 0 if you want 
%                                   to perform it manually (using 
%                                   mrAnatAverageAcpcNifti.m).
%  4) brakeaftert1            - [no longer in use]
%  5) brakeaftervisualization - [no longer in use]
%  6) channels                - [no longer in use]
%  7) check                   - Enter 1 if an interactive image is used to
%                                    detect movement, or 0 if not (default)
%  8) clobber                 - Enter 1 to overwrite existing data and
%                                    reprocess, or 0 (default) not to 
%                                    reprocess. As a mrQ_Set input it
%                                    determines whether to redo the whole
%                                    mrQ, but clobber also appears in
%                                    functions within mrQ, to detect
%                                    whether to reprocess a particular
%                                    function or sub-function.
%  9) coilweights
% 10) complexflag
% 11) fieldstrength          - Enter the field strength of the scanner.
% 12) freesurfer             - Enter 1 if FreeSurfer is available for use,
%                                     or 0 if not (default).
% 13) inputdata_seir         - Enter a string naming the directory of where 
%                                    the SEIR data is located
% 14) inputdata_spgr         - Enter a string naming the directory of where 
%                                    the SPGR data is located
% 15) interp                 - Enter 1 (default) to use the trilinear
%                                    method of interpolation, or 7 to use
%                                    the b-spline method. This is used in
%                                    aligning the SPGR images.
% 16) lsq                    - Enter 1 to run the least squares model for 
%                                    the T1-M0 fit, or 0 not to (default).
% 17) lw                     - Enter 1 (default) to run the linear weighted
%                                    model for the T1-M0 fit, or 0 not to.
% 18) mmpervox               - Enter a 3x1 vector of the image resolution.
%                                    The default is the NIfTI's resolution.
% 19) muti_coil_im
% 20) pdfit_method           - Enter 1 for T1 single coil (default), 2 for 
%                                    multi-coil, or 3 for both. This is 
%                                    used in the M0-to-PD fit.
% 21) permutation
% 22) polydeg                - Enter a scalar for the polynomial degree to 
%                                     fit to the Gain (default is 3)
% 23) proclus                - [no longer in use]
% 24) reset_raw_dir
% 25) rawdir                 - Enter a string naming the directory of the 
%                                     raw image, as defined in mrQ_Create
% 26) refim                  - Enter the path to a NIfTI file to use as a
%                                     reference image. If empty, the
%                                     default is that the SPGR with a
%                                     similar contrast to the T1-weighted
%                                     image will be selected, and the user 
%                                     will be asked to mark the ac-pc using 
%                                     mrAnatAverageAcpcNifti.m.
% 27) runfreesurfer
% 28) seir                   - Enter a cell of strings (e.g. {'001' 002'}) 
%                                     as the list of SEIR scan numbers. 
% 29) siemens                - [no longer in use]
% 30) siemensdicom           - [no longer in use]
% 31) skip                   - Enter a 1xN vector of scans to skip. Default
%                                   is not to skip any.
% 32) spgr                   - Enter a cell of strings (e.g. {'003' 004'}) 
%                                     as the list of SPGR scan numbers.  
% 33) sub                    - Enter a string as the name of subject.
% 34) sungrid                - Enter 1 to use the SunGrid calls, or 0 not
%                                     to (default).
% 35) useabs                 - [no longer in use]
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
        mrQ.autoacpc=0; 
        % if a ref image is given, the alignmant is mannual acpc
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
        
        % are there thresholding parameters for the PD build?
     case{'RepErrThreshold'}
          mrQ.RepErrThreshold=varargin;
     case{'PrcCutOff'}
         mrQ.PrcCutOff=varagin;
     case{'ErrorThresh'}
         mrQ.ErrorThresh=varargin;
    
    % threshold for spgr epi registration
    case{'AntsThresh', 'antsthresh'}
         mrQ.AntSQantTresh=varargin;    
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
