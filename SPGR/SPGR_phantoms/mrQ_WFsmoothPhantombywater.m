function [WFfile1,Gainfile1] =mrQ_WFsmoothPhantombywater(outDir,degree,CSFfile,pdfile,M0LFfile,BMfile,FSfile,SaveDir)
%
% mrQ_WFsmooth(outDir,CSFfile,pdfile,M0LFfile,BMfile,FSfile,SaveDir)
%
% # % the combined  PD from the fitted image box is not full PD map.
% becouse some image area have   no fits for different reason.
%the scombined and saved PD are therefore still with holls
% also some time the edged between the combined boxes can still exsist (in case of not a
% perfect scale estimation or boxs fits).
% the last step in done hear
% we  derive the gain by devide M0/pd=Gain for each coil. 
%we will asume that the gain must be smooth so we will keep it smooth all
%over the image
% and then get a a full PD by devide M0/gain=PD we will do it for each coil
% and then combine them again
%last we we will calculate the WF maps by normalized the PD to the CSF pd
%values

% INPUTS: 
% outDir         - The output directory - also reading file from there
%
% pdfile         - the combined pd image estimate by the boxes (using
%                    mrQ_BuildCoilsfitsPD)
%
% M0cfile        - The combined (lip angles)/aligned M0 data multi coil 4D
%
% BMfile         - The path to the brain mask
%
% FSfile         - The path to the free surfare segmentaiton mask (kliped by T1) 
%
% SaveDir        - the path for a saving directory
%
% OUTPUT         - the saved WF_map anf Gain file will be saved as a nii.gz file we
%                    output the path to those files.
%


%% CHECK INPUTS AND SET DEFAULTS
% load input files
if (~exist('degree','var')|| isempty(degree)),
    degree=3;
end;

if(~exist('SaveDir','var') ||  isempty(SaveDir) )
SaveDir=outDir;
end;

if(exist('CSFfile','var') &&  ~isempty(CSFfile) )
    disp(['Loading CSF data from ' CSFfile '...']);
    CSF = readFileNifti(CSFfile);
    CSF=double(CSF.data);

else

    CSFfile= fullfile(outDir,['csf_seg_T1.nii.gz']);
    disp(['trying to load CSF from ' CSFfile '...']);
    if(exist(CSFfile,'file'))
        disp(['Loading CSF data from ' CSFfile '...']);
        CSF = readFileNifti(CSFfile);
        CSF=double(CSF.data);
    else

        disp(['error , can not find the file: '  CSFfile]);
        error
    end
end;
if(exist('pdfile','var') &&  ~isempty(pdfile) )
    disp(['Loading pd data from ' pdfile '...']);
    PD = readFileNifti(pdfile);
    PD=double(PD.data);

else

    pdfile=fullfile(outDir,['M0_map_lsq.nii.gz']);
    disp(['trying to load PD from ' pdfile '...']);
    if(exist(pdfile,'file'))
        disp(['Loading PD data from ' pdfile '...']);
        PD = readFileNifti(pdfile);
        PD=double(PD.data);
    else

        disp(['error , can not find the file: '  pdfile]);
        error
    end
end;


% 
% if(exist('HMfile','var') &&  ~isempty(HMfile) )
%     disp(['Loading pd data from ' HMfile '...']);
%      HM = readFileNifti(HMfile);
%         HM=logical(HM.data);
% 
% else
% HMfile=fullfile(outDir,'HeadMask.nii.gz');
%     disp(['trying to load M0 from ' HMfile '...']);
%     if(exist(HMfile,'file'))
%         disp(['Loading  data from ' HMfile '...']);
%        HM = readFileNifti(HMfile);
%         HM=logical(HM.data);
%     else
% 
%         disp(['error , can not find the file: '  HMfile]);
%         error
%     end
% end;


if(exist('BMfile','var') &&  ~isempty(BMfile) )
    disp(['Loading pd data from ' BMfile '...']);
     BM = readFileNifti(BMfile);
     xform=BM.qto_xyz;
     mmPerVox=  BM.pixdim;
     BM=logical(BM.data);

else
BMfile=fullfile(outDir,'brainMask.nii.gz');
    disp(['trying to load M0 from ' BMfile '...']);
    if(exist(BMfile,'file'))
        disp(['Loading  data from ' BMfile '...']);
       BM = readFileNifti(BMfile);
        xform=BM.qto_xyz;
     mmPerVox=  BM.pixdim;
     BM=logical(BM.data);
    else

        disp(['error , can not find the file: '  BMfile]);
        error
    end
end;





%% normalize by CSF and save
if(exist('FSfile','var') &&  ~isempty(FSfile) )
    disp(['Loading segmentation data from ' FSfile '...']);
    fs = readFileNifti(FSfile);
    fs=fs.data;

else

FSfile=fullfile(outDir,['T1w_tissue.nii.gz']);
    disp(['trying to load  from ' FSfile '...']);
    if(exist(FSfile,'file'))
    disp(['Loading segmentation data from ' FSfile '...']);
          fs = readFileNifti(FSfile);
          fs=fs.data;

    else
        disp(['error , can not find the file: '  FSfile]);
        error
    end
end;

Imsz1=size(PD);

[Poly1,str] = constructpolynomialmatrix3d(Imsz1,find(ones(Imsz1)),degree);


%fit the polynials coefitents to the smooth B1 map
[params,gains,rs] = fit3dpolynomialmodel(PD,(fs==1),degree);

%reshape from a vector to 3D map
Gain = reshape(Poly1*params(:),Imsz1);

pd_corect=PD./Gain;



 wmV=mean(pd_corect(fs>1 & pd_corect>0));
  
      %THe ventrical CSF ROI is where the free surfare defind, where T1 is reisnable (see mrQ_CSF)  but it can't have a White matter  values and not double it (maybe worng for kids)
   CSF1= fs==1 & pd_corect>wmV & pd_corect< wmV*2;
      CSF1=CSF & pd_corect>wmV & pd_corect< wmV*2; 

 % CalibrationVal=median(pd_corect(find(CSF1)));
  % pd_corect=pd_corect./CalibrationVal;
% %  
%let find the scalre from the maxsimum of the csf values histogram
    [d dd]= ksdensity(pd_corect(find(CSF1)), [min(pd_corect(find(CSF1))):0.001:max(pd_corect(find(CSF1)))] );
% %  CalibrationVal=mean(PD(find(CSF)));
 CalibrationVal= dd(find(d==max(d)));% median(PD(find(CSF)));
 
 %lets calibrate the pd by the pd of the csf roi
    pd_corect=pd_corect./CalibrationVal;
    
% let cut outlayers
pd_corect(pd_corect<0)=0;
pd_corect(pd_corect>2)=2;
 %let save it
   WFfile1=fullfile(SaveDir,['WF_map.nii.gz']);
   dtiWriteNiftiWrapper(single(pd_corect), xform, WFfile1);
 Gainfile1=fullfile(SaveDir,['Gain.nii.gz']);
   dtiWriteNiftiWrapper(single(Gain), xform, Gainfile1);