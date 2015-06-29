function mrQ = mrQ_Set(mrQ,param,varargin,saveflag)
% Setting data in the mrQ structure for the mrQ anlysis
%
% The mrQ structure stores all the MRQ computations.  mrQ_Set is the main
% function to add values, images, computations etc. to this structure.  As
% inputs mrQ_Set takes in an mrQ structure, a parameter defined as a
% 'string' which is to be set, and a number of aditional arguments and
% values that are specific to that parameter.  The arguments for each
% parameter are described below.  mrQ_Set always returns the modified mrQ
% structure and save it unless saveflag=0.  The structure is created with defult values with mrQ_Create
% or in the first call.

%See also:  mrQ_Create

% Parameter list and associated arguments:
% saveflag          - the defult is yes; for no save make it 0
% param             - varargin the key parameters are
%   rawdir            - directory of raw image   string 'path' this define in mrQ_Create
%   seir              - list of SEIR scan number cell of string {'001' 002'}
%   SPGR              - list of SPGR scan number  cell of string{'001' 002'}
%   sub               - a name of subject string 'name'
%   clobber           - true or faluse to redo the anlysis
%   refim             - a path to an imge to use as arefernce string 'path'
%   check            - if a interactive image is used to ditect movment number - 1
%   proclus          -use the proclass number - 1
% sungrid          - use the sgesngrid calls -1
%  polydeg           - the polynomyals degree to fit to the Gain (defult 3)
%  'brakeaftervisualization' or 'viewbrake' - use that to start the data
%  visual the images wirh the check flag and then stop the code from
% To reverse sewpesipic steps in the fir
% seir_done                                1 to skip the SEIR t1 fit or 0 to redo it
% spgr_init_done                      1 to skip the spgr init stepor 0 to redo it
% spgr_coilweight_done           1 to skip wighting 0 to redo it
% spgr_t1fit_done                      1 to skip T1 fit SPGR 0 to redo it
% segmentation                            1 to skip segmentation 0 to redo it
%  calm0_done                           1 to skip M0 claculation 0 to redo it
% spgr_pdfit_done                     1 to skip PD fit0 to redo it
%spgr_pdbuild_done                1 to skip PD fit0 to redo it
%
% ...                  - check code for others
%
% example:
% mrQ = mrQ_Set(mrQ,,'raw','/biac2/wandell2/data/WMDevo/ambliopia/sub7/QuantativeImaging/20121102_3488')
% mrQ = mrQ_Set(mrQ,'sub','sub7');
% mrQ = mrQ_Set(mrQ,'SPGR',{ '0009' '0010' '0011' '0012' });
% mrQ = mrQ_Set(mrQ,'seir',{ '0004' '0005' '0006' '0007' });
% mrQ = mrQ_Set(mrQ,'proclass',1);


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
        
        error('Uknown mrQ parameter');
end


%save mrQ
if notDefined('saveflag')
    saveflag=1;
end

if saveflag==0
else
    save(mrQ.name,'mrQ');
end
