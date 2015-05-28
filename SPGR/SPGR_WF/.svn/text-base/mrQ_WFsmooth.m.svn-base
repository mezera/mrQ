function [WFfile1,Gainfile1] =mrQ_WFsmooth(outDir,CSFfile,pdfile,M0LFfile,BMfile,FSfile,SaveDir,Segfile)
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

    pdfile=fullfile(outDir,['PD_fitGboxMedian.nii.gz']);
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


%% load M0 image
if(exist('M0Lfile','var') &&  ~isempty(M0Lfile) )
    disp(['Loading pd data from ' M0Lfile '...']);
     M0 = readFileNifti(M0LFfile);
        M0=double(M0.data);

else
M0LFfile= fullfile(outDir,'/AligncombineCoilsM0.nii.gz');
    disp(['trying to load M0 from ' M0LFfile '...']);
    if(exist(M0LFfile,'file'))
        disp(['Loading M0 data from ' M0LFfile '...']);
        M0 = readFileNifti(M0LFfile);
        M0=double(M0.data);
 
    else

        disp(['error , can not find the file: '  M0LFfile]);
        error
    end
end;
%%
% Handle the Freesurfer Segmentation

% Load the segmentation nifti file - output from freesurfer
if(~exist('Segfile','var') ||  isempty(Segfile) )

segfile=fullfile(outDir,'t1_bet_seg.nii.gz');
end


fs = readFileNifti(segfile);

% Dimenstionality check 
if fs.dim(1)==size(BM,1) && fs.dim(2)==size(BM,2) && fs.dim(3)==size(BM,3)
else
    bb = mrAnatXformCoords(xform,[1 1 1;size(BM)]);
    fs.data = mrAnatResliceSpm(double(fs.data),inv(fs.qto_xyz),bb,mmPerVox,1);
end

fs = double(fs.data);

% Do some more dimensionality checks between the T1 and the freesurfer
% segmentation
if size(fs,1)==size(BM,1)+1;
    disp('freesurfer is different in size from the T1 we will clip the extra x voxels and hope its right')
    fs1 = fs; clear fs
    fs(1:size(T1,1),:,:)=fs1(1:size(T1,1),:,:); clear fs1;
end

if size(fs,2)==size(BM,2)+1;
    disp('freesurfer is differe in size from the data we clip the extra y voxels. we hope it right')
    fs1=fs;clear fs
    fs(:,1:size(BM,2),:)=fs1(:,1:size(BM,2),:); clear fs1;
end

if size(fs,3)==size(BM,3)+1;
    disp('freesurfer is differe in size from the data we clip the extra z voxels. we hope it right')
    fs1=fs;clear fs
    fs(:,:,1:size(BM,3))=fs1(:,:,1:size(BM,3)); clear fs1;
end

if size(fs,1)~=size(BM,1) || size(fs,2)~=size(BM,2) || size(fs,3)~=size(BM,3)
    error('The freesurfer segmentation file does not match the T1 data size!')
end



%% get teh estimate gain for each coil
mask=BM & PD>0 & fs==3 %(fs==2 | fs==41); %we will fit the smooth data when we have a good data (the white matter ).
  [XI YI]=meshgrid(1:size(PD,1),1:size(PD,2));
ZZI=zeros(size(PD));
V1=zeros(size(PD));
M01=zeros(size(M0));
%let run over coils and fit the Gain for each
for ii=1:size(M0,4)
    ii
   %the estimate gain for the coil
    V1=PD./M0(:,:,:,ii);
V1(~mask)=nan; 
V1(isinf(V1))=nan;

%intilaize the smooth gain and pd
ZZI=zeros(size(PD));
pd_corect=ZZI;
  
%loop over slices
for  i=1:size(PD,3)

        tmp=V1(:,:,i);
        
                %check that there is data in the slice
        wh=find(tmp>0);
        if  length(find(tmp>0))>100;
%find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
% estimate a smooth vertion of the data in the slice
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
             ZI = griddata(xg,yg,zg,XI,YI);         
            % put the result gain in the 3D gain image 
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            ZZI(:,:,i)=ZI;
            %put the corected pd in the 3D map
             pd_corect(:,:,i)=M0(:,:,i,ii).*ZI;
            clear ZI pdi
        end;

    end;
    %save the pd result for this coil
    M01(:,:,:,ii)=pd_corect;
   
end

%% avrage the coils
 % let make waites by SNR for each coil in each voxel
 
 W=sum(M0,4);
for i=1:size(M0,4)
    W1(:,:,:,i)=W;
end

W1=M0./W1;

[~, In]=sort(W1,4,'descend');


coilMask=nan(size(W1));
%let use the first 16 coils (so we wont avrage noise). if there are less
%then we use all the coil we got
for i=1:min(size(In,4),16)
   T=In(:,:,:,i); 
    for ii=1:min(size(In,4)) %
   TT=find(T==ii);
   TTT= coilMask(:,:,:,ii)    ;
   TTT(TT)=1;
   coilMask(:,:,:,ii)=TTT;
        
    end
end
clear In T TT TTT
% let keep only the value we like (high SNR) the other will be nans 
II=M01.*coilMask;

II(find(II<0))=nan; % if we have negativ value it fits of noise...

pd_corect=nanmedian(II,4);
clear II coilMask;

%W=W1; 
%we combine the PD of each coil according to the waits
%  pd_corect=median(M01,4);
    %pd_corect1=(sum(M01.*W1,4));
   % the median is defently better and we don't get the remained coils biases but we loss some SNR. maybe it is good to treshold the coil by there precental signal below some value is better nit to use
 %  and then the SNR will be better (?)
   
    M0_=sqrt(sum(M0.^2,4));
     clear M0




%% normalize by CSF and save
if(exist('FSfile','var') &&  ~isempty(FSfile) )
    disp(['Loading segmentation data from ' FSfile '...']);
    fs = readFileNifti(FSfile);
    fs=fs.data;

else
FSfile = fullfile(outDir,'T1w_tissue.nii.gz')
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

% find wm mask for pd
   wmV=mean(pd_corect(fs==2 & pd_corect>0));
  
      %THe ventrical CSF ROI is where the free surfare defind, where T1 is reisnable (see mrQ_CSF)  but it can't have a White matter  values and not double it (maybe worng for kids)
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
    Gain=M0_./pd_corect;
% let cut outlayers
pd_corect(pd_corect<0)=0;
   pd_corect(pd_corect>2)=2;

   %let save it
   WFfile1=fullfile(SaveDir,['WF_map.nii.gz']);
   dtiWriteNiftiWrapper(single(pd_corect), xform, WFfile1);
 Gainfile1=fullfile(SaveDir,['Gain.nii.gz']);
   dtiWriteNiftiWrapper(single(Gain), xform, Gainfile1);

  return 
%%  old code 


% 
%  wmV=mean(pd_corect(fs==2 & pd_corect>0));
%   
% 
% 
% %% fit smooth bias feild
% 
%    mask=BM & PD>0;
%     V1=PD./M0;
% V1(~mask)=nan; 
% V1(isinf(V1))=nan;
% 
%  ZZI=zeros(size(PD));
% pd_corect=ZZI;
%     [XI YI]=meshgrid(1:size(PD,1),1:size(PD,2));
% 
% 
% for  i=1:size(PD,3)
% 
%         tmp=V1(:,:,i);
%         wh=find(tmp>0);
% 
%         if  length(find(tmp>0))>100;
% 
%             [x,y] = ind2sub(size(tmp),wh);
%             z=double(tmp(wh));
% 
%             [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',20);
%              ZI = griddata(xg,yg,zg,XI,YI);         
%             ZI=rot90(ZI);
%             ZI = flipdim(ZI,1);
%             ZZI(:,:,i)=ZI;
%              pd_corect(:,:,i)=M0(:,:,i).*ZI;
%             clear ZI pdi
%         end;
% 
%     end;
% 
% 
% 
% mask=BM & PD>0;
% 
% M01=zeros(size(M0));
% for i=1:size(M0,4)
%     i
%     tic
%     V1=PD./M0(:,:,:,i);
% V1(~mask)=nan; 
% V1(isinf(V1))=nan;
% 
%  ZZI=zeros(size(PD));
% pd_corect=ZZI;
%     [XI YI]=meshgrid(1:size(PD,1),1:size(PD,2));
% 
% 
% for  i=1:size(PD,3)
% 
%         tmp=V1(:,:,i);
%         wh=find(tmp>0);
% 
%         if  length(find(tmp>0))>100;
% 
%             [x,y] = ind2sub(size(tmp),wh);
%             z=double(tmp(wh));
% 
%             [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',20);
%              ZI = griddata(xg,yg,zg,XI,YI);         
%             ZI=rot90(ZI);
%             ZI = flipdim(ZI,1);
%             ZZI(:,:,i)=ZI;
%              pd_corect(:,:,i)=M0(:,:,:,i).*ZI;
%             clear ZI pdi
%         end;
% 
%     end;
%     
%     M01(:,:,:,i)=pd_corect;
%     toc
% end
% 
%      W=sum(M0,4);
%  
% for i=1:size(M0,4)
%     W1(:,:,:,i)=W;
% end
% 
% W1=M0./W1;
%  
% 
%     M02=(sum(M01.*W1,4));
%       clear M0 W1 M01
%          M0=M02;
% clear M02   
% 
% 
% 
