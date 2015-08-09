function [mrQ] = mrQ_WF(mrQ,dataDir,T1file,PDfile,Gainfile,M0file);

    
    %% II. Load aligned data
if notDefined('dataDir');
    dataDir = mrQ.spgr_initDir;
end
outFile  = fullfile(dataDir,'dat_aligned.mat'); %without coilWeights data

disp(['Loading aligned data from ' outFile '...']);

load(outFile);

if notDefined('T1file')
       [ T1file,~,~]=mrQ_get_T1M0_files(mrQ,1,0,0);
end
T1=readFileNifti(T1file);
mmPerVox = T1.pixdim;
T1=T1.data;

if notDefined('PDfile')
    load(mrQ.opt_logname);
    PDfile=opt.PDfile;
end

% need to check which M0 map will be here...

if notDefined('M0file')
    M0file=opt.M0file_Org;
end
 M0=readFileNifti(M0file);
 outMm=M0.pixdim(1:3);
 M0=M0.data; 
 MSZ=size(M0);
%     %M0=sqrt(sum(M0.^2,4));
%     M0=(sum(M0,4));

WFmask=ones(size(T1));
%% I
% find area of high T1=4.3 +-0.1 around the vetricles

% clip around the ventricals:
if notDefined('boxsize')
    boxsize(1)=30;
    boxsize(2)=40;
    boxsize(3)=20;
end
sz=size(WFmask); szH=round(sz./2);
XX=boxsize(1)./round(mmPerVox(1));
YY=boxsize(2)./round(mmPerVox(2));
ZZ=boxsize(3)./round(mmPerVox(3));

WFmask(szH(1)+XX:end,:,:)=0;
WFmask(1:szH(1)-XX,:,:)=0;

WFmask(:,1:szH(2)-YY,:)=0;
WFmask(:,szH(2)+YY:end,:,:)=0;

WFmask(:,:,1:szH(3)-ZZ)=0;
WFmask(:,:,szH(3)+ZZ:end)=0;

% find areas within the ventricles area with high value of T1
WFmask= WFmask & T1<=4.4 & T1>=4.2 ;


%% II
%loop over raw (aligned and combined) data and find the flip angle(s) with
%the best SNR. it is hard to measure noise. so we will use the signal in WM
%as a reference because it's a place with good signal across flip angles)
fa=SPGR_niiFile_FA;
for ii=1:length(s.imdat)
    
    raw=readFileNifti(mrQ.AlignedSPGR{ii};
    rawDat=raw.data;
    
    %     WMmask=? 
    %      per file??

    WMrawDat=rawDat(WMmask);
    WMdata(ii)=mean(WMrawDat(:));
    
end

FAloc=find(WMdata==max(WMdata));
bestFA=fa(FAloc);

%% III
% exclude M0 outliers in the ROI

WFmask=WFmask & M0<prctile(M0(WFmask),99);


%%
%4. calculate  PD within the selected ROI and flip angle(s) (using the SPGR
%equation, T1, and gain [note gain for multi coil is define for the sum of
%multi coils, and for the combine as arrived from the scanner])




%%
%5 scale this roi to have pd of 1. --> this is our global scale

%%
%6. aplly to PD images to make WFfile --> done
