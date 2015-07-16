function mrQ=mrQ_segmentT1w2tissue(mrQ,BMfile,T1file,t1wfile,outDir,csffile,boxsize)
% mrQ=mrQ_segmentT1w2tissue(mrQ,[],mrQ.SegInfo.T1wSynthesis_T1);
% 
% The function segment by FSL to three tissue take the CSF tissue restrict
% it by the T1 values the CSF is also restricted to be in the center of the
% brain in a box of ~ 60x80x40mm asuming the brain is in ACPC space this is
% where the ventrical are.
% 
% 
if~exist('outDir')
outDir = mrQ.spgr_initDir;
end

if ~exist('T1file','var') || isempty(T1file)
   [T1file, ~,~]=mrQ_get_T1M0_files(mrQ,1,0,0);
end

disp(['Loading T1 data from ' T1file '...']);
T1 = readFileNifti(T1file);
xform = T1.qto_xyz;
mmPerVox = T1.pixdim;
T1 = double(T1.data);


%     T1file = fullfile(outDir,'T1_map_lsq.nii.gz');
%     disp(['trying  to load T1 from ' T1file '...']);
%     if(exist(T1file,'file'))
%         disp(['Loading T1 data from ' T1file '...']);
%         T1 = readFileNifti(T1file);
%         T1 = double(T1.data);
%     else
%         T1file = mrvSelectFile('r','*.nii.gz','Select T1 File');
%         % disp(['error , can not find the file: '  T1file]);
%         %error
%     end
%end

if notDefined('BMfile')
    BMfile = fullfile(outDir,'brainMask.nii.gz');
    if ~exist(BMfile,'file')
        BMfile = mrvSelectFile('r','Select the Brain Mask');
    end
end
BM = readFileNifti(BMfile); BM=BM.data;


if notDefined('t1wfile')
    t1wfile=mrQ.SegInfo.T1wSynthesis_T1;
    if ~exist(t1wfile,'file')
        t1wfile = mrvSelectFile('r','Select the T1 weighted image');
    end
end




%% Load Data


% Load the Analysis info file and append it with the info provided
% infofile = fullfile(outDir,'AnalysisInfo.mat');
% load(infofile);

% AnalysisInfo.GainPolydegrees = degrees;
% AnalysisInfo.M0forWF = M0file;


% Load the brain mask from the outDir
% disp(['Loading brain Mask data from ' BMfile '...']);
% brainMask = readFileNifti(BMfile);
% xform = brainMask.qto_xyz;
% mmPerVox = brainMask.pixdim;
% brainMask = logical(brainMask.data);


%% FSL
% % % run bet of fsl
% % % fprintf('\n running fsl bet ...\n');
% % % 
  t1w1file=mrQ.SegInfo.T1wSynthesis_MOT1;
segfile=fullfile(outDir,'t1_bet.nii.gz');
% eval(['! bet ' t1wfile ' ' segfile ' -m -f .4'])
% % % 
eval(['! bet ' t1w1file ' ' segfile ' -m -f .4'])


fprintf('\n fsl fast ...\n');
% eval(['! fast '  TforSeg_file ])
eval(['! fast '  segfile ])

 segfile=fullfile(outDir,'t1_bet_seg.nii.gz');
%  segfile=fullfile(outDir,'TforSeg_seg.nii.gz');


% AnalysisInfo.T1forWF        = T1file;
% AnalysisInfo.WFdate         = date;
% AnalysisInfo.T1wSeg         = segfile
% save(infofile,'AnalysisInfo');


fprintf('\n T1 map and ventrical location restrictions... \n');

% T1 cliping of csf
seg=readFileNifti(segfile);

%check if fsl conventions correspond predictions regarding segmented T1
%  seg(GM)=1, seg(WM)=2, seg(CSF)=3
tmpseg=seg.data; tmpseg(~logical(BM))=0;
clus1=mode(T1(tmpseg==1)); clus2=mode(T1(tmpseg==2)); clus3=mode(T1(tmpseg==3));
csfClus=find([clus1,clus2,clus3]==max([clus1,clus2,clus3]));
WMClus=find([clus1,clus2,clus3]==min([clus1,clus2,clus3]));
GMClus=find([clus1,clus2,clus3]==median([clus1,clus2,clus3]));
clear tmpseg;
 if ~(csfClus==3 & GMClus==1 & WMClus==2) 
     warn=sprintf('The fsl segmentation is possibly false, possibly due to extreme values of T1.\n It is recommended to segment te CSF elsewhere and insert the csf nifti as input.\n Segmentation can be done either mannualy or using freesurfer.\n ');
     error(warn)
 end

CSF1=zeros(size(T1));CSF1(seg.data==3)=1;

%cliping the center box

if notDefined('boxsize')
    boxsize(1)=30;
    boxsize(2)=40;
    boxsize(3)=20;
end
sz=size(CSF1); szH=round(sz./2);
XX=boxsize(1)./round(mmPerVox(1));
YY=boxsize(2)./round(mmPerVox(2));
ZZ=boxsize(3)./round(mmPerVox(3));

CSF1(szH(1)+XX:end,:,:)=0;
CSF1(1:szH(1)-XX,:,:)=0;

CSF1(:,1:szH(2)-YY,:)=0;
CSF1(:,szH(2)+YY:end,:,:)=0;

CSF1(:,:,1:szH(3)-ZZ)=0;
CSF1(:,:,szH(3)+ZZ:end)=0;

if notDefined('csffile')
    csffile = fullfile(outDir, 'csf_seg_T1.nii.gz');
end

CSFtmp=CSF1;

CSF1= CSF1 & T1>4. & T1< 5;
dtiWriteNiftiWrapper(single(CSF1), xform, csffile);

if length(find(CSF1))<200
           fprintf(['\n Warning: We could find only ' num2str(length(find(CSF1))) ' csf voxel this make the CSF WF estimation very noise cosider to edit csf_seg_T1.nii.gz fiel see below \n']);
end
% we save also ROI that  will include more voxel in the CSF ROI but will be
% more permisive to patial voulume. this should be use if CSF is hard hard
% to estime form small number of voxel. this can happen with subject with
% small vertricals and low resulotion scans.
CSFtmp= CSFtmp & T1>3;

%length(find(CSF1));
csffile1 = fullfile(outDir, 'csf_seg_T1_large.nii.gz');
dtiWriteNiftiWrapper(single(CSFtmp), xform, csffile1);
clear CSFtmp
%% seg(GM)=1, seg(WM)=2, seg(CSF)=3
% mask = zeros(size(brainMask));
mask = zeros(size(T1));
mask = double(mask);
mask(find(CSF1)) = 1;

% Set up the white matter mask fsl
% wm = zeros(size(brainMask));
wm = zeros(size(T1));
wm(seg.data==2)  = 1;

wm         = logical(wm);
[d dd] = ksdensity(T1(wm),(min(T1(wm)):0.01:max(T1(wm))));

M = dd(find(d==max(d)));
% M = mean(T1(wm)) ;
% S = std(T1(wm));

wm = wm &  T1>(M-0.03) & T1<(M+0.03);
mask(find(wm)) = 2;

% Set up the gray matter mask fsl
% cortex = zeros(size(brainMask));
cortex = zeros(size(T1));
cortex(seg.data==1) =1;
cortex         = logical(cortex);

[d dd] = ksdensity(T1(cortex), (min(T1(cortex)):0.01:max(T1(cortex))) );
M      = dd(find(d==max(d)));
cortex = cortex & T1 > (M-0.03) & T1 < (M+0.03);

mask(find(cortex)) = 3;
filefsl = fullfile(outDir,'T1w_tissue.nii.gz');
dtiWriteNiftiWrapper(single(mask), xform, filefsl);

mrQ.csf_large=csffile1;
mrQ.csf=csffile;
mrQ.T1w_tissue=filefsl;





%% run ANTS segmentation using the brainmask aand the Atropos command
% %  segfile=fullfile(outDir,'t1_seg.nii.gz');
% %  seg_prob_file=fullfile(outDir,'t1_seg_prob.nii.gz');

% % cmAtropos=['/home/ANTs-2.1.0-Linux/bin/Atropos -d 3 -a ' t1wfile ' -o [' segfile ' , ' seg_prob_file '] '...
% % '-c [10,0.001] -m [0.1,1x1x1] -i kmeans[4] -x ' BMfile ];
% % 
% % % Run the command in unix and get back status and results: 
% %     [status, result] = system(cmAtropos); 
%%or instead of BET USE THE BRAIN maSK
% TforSeg_file=fullfile(outDir,'TforSeg.nii.gz');
% TforSeg=readFileNifti(t1wfile); TforSeg=TforSeg.data;
% TforSeg(BM==0)=0;
% dtiWriteNiftiWrapper(single(TforSeg),xform,TforSeg_file);