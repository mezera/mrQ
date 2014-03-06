function mrQ_run(mrQfileName,clobber)
%mrQ_run(RawDir,arrangeRawFlag,SEIR_seriesNumbers,SPGR_seriesNumbers,refIm,sub,freesurfer,channels,useNiftiFlag,alignFlag,complexFlag,useAbs,mmPerVox,interp,skip,coilWeights,clobber)
% mrQ_run(RawDir,arrangeRawFlag,SEIR_seriesNumbers,SPGR_seriesNumbers,refIm,sub,freesurfer,channels,useNiftiFlag,alignFlag,complexFlag,useAbs,mmPerVox,interp,skip,coilWeights,clobber)
%
% # a batch call to run all the mrQ fits in one click
% # the batch  run over mrQ functions to claculate the maps (T1,PD TV WF
% SIR) it will also run Ants SGE and freesurfer softwhere. all of those and
%%the mrVista reposotorry must be set and updated
%
% INPUT VARIABLES of the structure mrQ:
%   RawDir            - The path of the directory containing the raw MR
%                       data.the function assume the data are arrange in the same way as it is in the Stanford CNI magnet.
%                       the data should be saved as niifti files and a dicom directory with dicom are also exist.
%                       This convention may be different for each scanner.
%                       It is good idea to cheack this and maybe arrange the data with this conventions or edit the first function that call the raw directory ( mrQ_arrangeData.m)
%
%   arangeRawFlag     - Needs to be 1 to move every file from the original
%                       magnet directory (dataDir) to a raw directory
%                       inside it (if this was done set = 0, or empty).
%
%   SEIR_seriesNumbers- A cell array containing the SEIR_epi series numbers
%                       from dataDir/raw.  (if this input exists a SEIR dir
%                       will be created with data and fitT1_GS directories
%                       inside it). The scan dicom will be extracted from
%                       the data Dir. If SEIR_epi already exists and a new
%                       one will be created.
%
%   SPGR_seriesNumbers- A cell array with the SPGR series numbers from
%                       dataDir/raw. This input creates a SPGR directory
%                       with data inside it. The scan dicom will be
%                       extracted from the data Dir. If SPGR already exists
%                       a new one will be created.
%
%  refImg              Different ref images can be used as an input (refImg is a
%                      path to a nifti image). If there is no refImg (refImg is
%                      empty) then the SPGR with a similar contrast to the T1
%                      weighted image will be selected and the user will be
%                      asked to mark the ac/pc using mrAnatAverageAcpcNifti
%
%       sub            - Subject name used for different file names
%                       like the SGE. if the name is not porvided we will take it for the
%                        raw directory
%
%   freesurfer         -an exssisting freesurfer segmentation aparc+aseg in
%                       nii.gz format. if this is not provied then freesurfer segmentation on a
%                       syntetic T1w image that we will calculate from the data. or on the refImg if it been provided. the freesurfer
%                       segmentation may take 24 hours. freesurfer segmentation may failed and
%                       the user need to cheack if it is.
%
%   channels          - A 1xN array (where N is the number of series). If
%                       any of the SPGR data has multi-channel dicoms make
%                       the series number equal to the number of channels.
%                       The default is zero.
%
%   useNiftiFlag      - A 1xN logical array (where N is the number of
%                       series). If it is set to 1 the code will use the
%                       nii.gz file and the channel variable input is not
%                       used. defult is ones
%
%   alignFlag          - if to align the SEIR data set the defult is 1 (yes).
%
%  complexFlag         - if the SEIR data is complex and not only magnitude the defult is 0 (no).
%
%   useAbs              -if the SEIR data is complex but only magnitude need to be used defult is 0 (no).
%
%   mmPerVox           - The resolution at which you want to resample the data.
%                        This is a 3X1 (1x3?) vector. If empty, the dicom
%                        resolution will be used-this does not have to be the
%                        native scan size as the magnet output % is zeroed. The
%                        saved directory will have the resolution in its name.
%
%  interp              -Interpolation method. [Default = 1]
%                       1 = trilinear,
%                       7 = b-spline (resampling algorithm)
%
%  skip                -you can skip any of the scans in spgrDir if you want by
%                       passing in a 1xn vector of scans to skip.
%
%  coilWeights         -determine optimal coil weights and save them out.
%                       Default = true; when we use more then 8 chanels
%                       coils
%
%  clobber              -Overwrite existing data and reprocess. [default = false]
%
%useNiftiFlag      - A 1xN logical array (where N is the number of
%                       series). If it is set to 1 the code will use the
%                       nii.gz file and the channel variable input is not
%                       used. *** WHY USE THIS OR DON'T? ***
%

%output % save a mrQparm file the document the anaysis. the different file
%make fiels and directory. the run end when a map directory with five map
%is crated including T1 WF TV SIR and VIP maps file
%
%Example :
% make or load exsist structure
%
%mrQ=mrQ_Create('/biac2/wandell2/data/WMDevo/ambliopia/sub7/QuantativeImaging/20121102_3488');
% set the must feileds
% mrQ=mrQ_Set(mrQ,'SEIR',{'0004' '0005' '0006' '0007'  })
% mrQ=mrQ_Set(mrQ,'SPGR',{'0009' '0010' '0011' '0012'  })
% mrQ=mrQ_Set(mrQ,'sub','amb_7')
%
%  runn or run the mrQ fit
% mrQ_run(mrQ)



%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~exist(mrQfileName,'file') )
    error(['cant find ' mrQfileName ])
else
    load (mrQfileName)

end



%% param for the mrQ_arrangeData

if ~isfield(mrQ,'RawDir')
    disp('Raw data need to be provided');
    error
end;


if ~isfield(mrQ,'SPGR_seriesNumbers')
    mrQ.SPGR_seriesNumbers=[];
end

% This will use the nii.gz file and not the dicom in the SPGR dir
% (should be the defult if they avilable)
if ~isfield(mrQ,'useNiftiFlag')
    mrQ.useNiftiFlag(1:length(mrQ.SPGR_seriesNumbers))=1;
end

sub = mrQ.sub;




%% arange data form the magnet

% most subject have this allready, we need to check before we run it should
% come form the user


if isfield(mrQ,'Arange_Date') && ~isempty(mrQ.Arange_Date)

    fprintf([ '\n  Arranged data found: ' mrQ.Arange_Date ' !  \n']);

else
    fprintf('\n  Arranging the data ...\n  ')
    [mrQ]=  mrQ_arrangeData(mrQ.RawDir,mrQ.arrangeRawFlag,mrQ.SEIR_seriesNumbers,mrQ.SPGR_seriesNumbers,mrQ.channels,mrQ.useNiftiFlag,mrQ);
    fprintf('arange data is done!\n');
    save(mrQ.name,'mrQ');
end

% usege of  SUnGrid is advised and is the defult.it can be turn off with
% mrQSet : mrQ=mrQ_Set(mrQ,'sungrid',0)

if ~isfield(mrQ,'SunGrid');
    mrQ.SunGrid = 1;
end


%% fit SEIR
% load(name); %we load the SEIR dir that was saved by mrQ_arrangeData before


if isfield(mrQ,'SEIR_done');
else
    mrQ.SEIR_done=0;
end

if (mrQ.SEIR_done==0);

    %keep track of the variable we use  for detail see inside the function
    [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag,mrQ.complexFlag,mrQ.useAbs);

    [mrQ]=mrQ_fitSEIR_T1(mrQ.SEIRepiDir,[],[],0,mrQ);
    mrQ.SEIR_done=1;
    save(mrQ.name,'mrQ');
    fprintf('fit SEIR  - done!');

else
    fprintf('\n  load fit SEIR data ! \n');

end

%% intiate and  Align SPGR
%  param for Align SPGR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load(name);

if isfield(mrQ,'SPGR_init_done');
else
    mrQ.SPGR_init_done=0;
end

if     mrQ.SPGR_init_done==0

    %keep track of the variable we use  for detail see inside the function
    [~, ~, ~,~,~, mrQ]=mrQ_initSPGR(mrQ.SPGR,mrQ.refIm,mrQ.mmPerVox,mrQ.interp,mrQ.skip,[],mrQ);
    mrQ.SPGR_init_done=1;

    save(mrQ.name,'mrQ');
    fprintf('\n  init SPGR - done!           \n');
else
    fprintf(' \n load init SPGR data            \n');

end


if isfield(mrQ,'brakeAfterVisualization');
else
    mrQ.brakeAfterVisualization=0;
end

if     mrQ.brakeAfterVisualization==1
    return
end

%% Save out a data structure with the optimal coil weighting applied



if isfield(mrQ,'SPGR_coilWeight_done');
else
    mrQ.SPGR_coilWeight_done=0;
end



if  (mrQ.coilWeights==1 && mrQ.coilNum(1)>8  && mrQ.SPGR_coilWeight_done==0)

    fprintf('\n Determining optimal coil weighting...\n');
    % Should this return the new structure with the weighting applied?
    [mrQ.AligndSPGR]=mrQ_multicoilWeighting(mrQ.spgr_initDir,mrQ.SPGR_niiFile,mrQ.SPGR_niiFile_FA,mrQ);
    mrQ.SPGR_coilWeight_done=1;
    fprintf('\n SPGR  coil weighting - done!               \n');

else

end
save(mrQ.name,'mrQ');


%%  Fit SPGR PD

if isfield(mrQ,'SPGR_T1fit_done');
else
    mrQ.SPGR_T1fit_done=0;
end

%vclover is implamented inside (we canadd rthis to the inputs
if (mrQ.SPGR_T1fit_done==0);
    [mrQ.AnalysisInfo]=mrQfit_T1M0_ver2(mrQ);
    mrQ.SPGR_T1fit_done=1;

    save(mrQ.name,'mrQ');
    fprintf('\n fit T1 SPGR  - done!              \n');
else
    fprintf('\n load fited  SPGR T1                \n');

end
save(mrQ.name,'mrQ');


%
% if we plan to run the freesurfer as part of the batch then a good time to
% do it is after mrQ_T1wSynthesis.m  is done (line 357 in mrQfit_T1M0_ver2). this function crate the T1
% wighted image that can be used as an input for freesurfer.

if isfield(mrQ,'brakeAfterT1');
else
    mrQ.brakeAfterT1=0;
end

if     mrQ.brakeAfterT1==1
    fprintf('\n brake T1 map  are done              \n');

    return
end

%% prefer to PD fit 1. get a segmentation (need freesurfer output) 2. get CSF; 3.make a M0 fies for the coils

%. Segmentaion and CSF
if isfield(mrQ,'segmentaion');
else
    mrQ.segmentaion=0;
end

if mrQ.segmentaion==0;


    % run Free surfare
    if (mrQ.runfreesurfer==1)
        mrQ=mrQ_Complitfreesurfer(mrQ);

        mrQ.segmentaion=1;
        % use an uploaded freesurafre nii.zg
    elseif isfield(mrQ,'freesurfer');
        [mrQ.AnalysisInfo]=mrQ_CSF(mrQ.spgr_initDir,mrQ.freesurfer,[],mrQ.AnalysisInfo);

        mrQ.segmentaion=1;

    else
        % Segment the T1w by FSL (step 1) and get the tissue mask (CSF WM GM) (step 2)
        mrQ=mrQ_segmentT1w2tissue(mrQ);
        mrQ.segmentaion=1;

    end
    save(mrQ.name,'mrQ');
end
%%
%3. coils MO
if isfield(mrQ,'calM0_done');
else
    mrQ.calM0_done=0;
end

if (mrQ.calM0_done==0);
    fprintf('\n calculate M0 for each coil               \n');

    %build a multi coil M0 image for the coils raw data and then fitted T1
    [mrQ.M0combineFile] = mrQ_multicoilM0(mrQ.spgr_initDir,[],[],mrQ.SPGR_niiFile,mrQ.SPGR_niiFile_FA,mrQ);

    mrQ.calM0_done=1;
    save(mrQ.name,'mrQ');
else
    fprintf('\n load M0 of each coil               \n');

end

%

%% fit PD

if isfield(mrQ,'SPGR_PDfit_done');
else
    mrQ.SPGR_PDfit_done=0;
end


if mrQ.SPGR_PDfit_done==0;
    fprintf('\n calculate PD from the M0 of all the coil               \n');


    %send the a call for the grid to fit  PD
%  [mrQ.opt]=mrQ_fitPD_multicoil(mrQ.spgr_initDir,1,[],mrQ.PolyDeg,[],mrQ.sub,mrQ.proclass);
 if isfield(mrQ,'opt_logname');
else
[mrQ.opt_logname]=mrQ_PD_multicoil_RgXv_GridCall(mrQ.spgr_initDir,mrQ.SunGrid,mrQ.proclus,mrQ.sub,mrQ.PolyDeg,[],[],[],[],[],[],[],[],[],mrQ);
save(mrQ.name,'mrQ');
mrQ_fitM0boxesCall(mrQ.opt_logname,mrQ.SunGrid,mrQ.proclus);


 end

    if mrQ.SunGrid==1;
        %check for SGE is done  before you move on (while loop)
        mrQ.SPGR_PDfit_done=mrQ_Gridcheack(mrQ.opt_logname,mrQ.SunGrid,mrQ.proclus);
    else
        mrQ.SPGR_PDfit_done=1;
    end

    save(mrQ.name,'mrQ');
else
    fprintf('\n load the claculted PD from the coils  M0               \n');

end



% bild the grid PD fits
if isfield(mrQ,'SPGR_PDBuild_done');
    %we need a check for the GE
else
    mrQ.SPGR_PDBuild_done=0;
end

if (mrQ.SPGR_PDBuild_done==0 && mrQ.SPGR_PDfit_done==1)
    fprintf('build the WF map   form PD fit            ');


    mrQ.opt=mrQ_buildPD(mrQ.opt_logname);

    %[mrQ.PDcombineInfo]=mrQ_BuildCoilsfitsPD(mrQ.spgr_initDir,mrQ.proclass);
    mrQ.SPGR_PDBuild_done=1;
    save(mrQ.name,'mrQ');
else
    fprintf('load the WF map               ');

end

% calculate VIP TV and SIR

if (mrQ.SPGR_PDBuild_done==1)
    fprintf('\n calculate VIP TV SIR form T1 and WF maps               \n');

    [mrQ.AnalysisInfo] = mrQ_VIP(mrQ);
    save(mrQ.name,'mrQ');
end

%% Organize the OutPut  directory
mrQ=mrQ_arangeOutPutDir(mrQ);

%done
mrQ.AnalysisDone=1;
mrQ.AnalysisDoneDate=date;
%save
save(mrQ.name,'mrQ');
fprintf('\n done !!! \n')
%%
