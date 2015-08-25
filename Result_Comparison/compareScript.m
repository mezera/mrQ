

T1wfile_1='T1wFilePath_FirstScan';
 T1wfile_2='T1wFilePath_SecondScan';
 outPutDir='Directory for Output';
 morefiles2_32ch={'PathForWFmap_1','PathForWFmap_2', 'PathForBrainMaskFile_1', 'PathForBrainMaskFile_2'};

% % register the files to one another: 
 WarpFiles=mrQ_ANTS_warp_SPGR2SPGR(T1wfile_1,T1wfile_2,outPutDir,morefiles2_32ch);

%% load the files: 
%T1
T1S1=readFileNifti(T1wfile_1);
T1S2=readFileNifti(fullfile(outPutDir,'Warp_T1w2_to_T1w1.nii.gz'));
% brain masks:
BMS1=readFileNifti(fullfile(outPutDir,WarpbrainMask1.nii.gz));
BMS2=readFileNifti(fullfile(outPutDir,WarpbrainMask2.nii.gz));
% WF files:
WF1=readFileNifti(fullfile(outPutDir,WarpWF1.nii.gz));
WF2=readFileNifti(fullfile(outPutDir,WarpWF2.nii.gz));

%% make a brainmask containing only brain tissue 

BM = logical(BMS2.data) & logical(BMS1.data);

% This region is cut off at the edge of scan box
BM(:,:,1:52)   = 0;
BM(:,:,144:end)= 0;

% We want T1 that is greater then 0.5 and smaller then 2.5. % this range
BM=BM & T1S1.data>0.7 & T1S1.data<2. & T1S2.data>0.7 & T1S2.data<2. ;


%% The coefficient of determination (R^2)
CV = (calccod(WF1.data(BM),WF2.data(BM))/;100)


%% Plot a  2D histogram of the error

% We plot the MTV = 1 - WF the macromolecular tisuue volume fraction.
nBins = 155;
[n,x,y] = mrQ_hist2d(WF1.data(BM),WF2.data(BM),nBins);

maxN = ceil(max(n(:))/10)*10;


mrvNewGraphWin;
image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256)));
identityLine(gca);
axis square xy; axis image;
%xlim([0.5 2]); ylim([0.5 2])
xlim([0. 0.7]); ylim([0. 0.7])

ylabel('WF  S1' ,'FontSize',16);
xlabel('WF  S2' ,'FontSize',16);
set(gca,'FontSize',16)
title( [' R^2= ' num2str(CV)])
grid on

showMontage(input2(:,:,90)); caxis([0 4]);axis off
showMontage(input1(:,:,90));caxis([0 4]);axis off


showMontage((input1(:,:,90)-input2(:,:,90))./input1(:,:,90));caxis([-0.5 0.5]);axis off





