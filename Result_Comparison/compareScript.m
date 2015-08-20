

% T1wfile18ch='/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1297_8ch_1mm3/output/OutPutFiles_1/T1w/T1w1.nii.gz';
% T1wfile2_32ch='/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1208_32ch_1mm/output/OutPutFiles_1/T1w/T1w1.nii.gz';
% outPutDir='/home/shai.berman/Documents/Code/AFQ_test/008_AM/WF_8ch_32ch_Comparisons';
% morefiles2_32ch={'/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1208_32ch_1mm/output/OutPutFiles_1/BrainMaps/WF_map.nii.gz', ... 
% '/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1208_32ch_1mm/output/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz'};
% 
% WarpFiles=mrQ_ANTS_warp_SPGR2SPGR(T1wfile18ch,T1wfile2_32ch,outPutDir,morefiles2_32ch);

%32 Chanel
% WF32chWarp=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/WF_8ch_32ch_Comparisons/WarpWF_map.nii.gz');
% BM32chWarp=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/WF_8ch_32ch_Comparisons/WarpbrainMask.nii.gz');

% WF32ch=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/WF_8ch_32ch_Comparisons/WF_32ch_4.2_4.7_kdensity.nii.gz');
% BM32ch=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1208_32ch_1mm/output/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz');
% T132=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1208_32ch_1mm/output/OutPutFiles_1/T1w/T1w1.nii.gz');
% T132.data=10.*T132.data;
%  %8 Chanel:
% WF8ch=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/WF_8ch_32ch_Comparisons/WF_8ch_4.2_4.7_kdensity.nii.gz');
% BM8ch=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1297_8ch_1mm3/output/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz');
% T18=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1297_8ch_1mm3/output/OutPutFiles_1/T1w/T1w1.nii.gz');
% T18.data=10.*T18.data;

WFS1=readFileNifti('/home/data/Stanford_Data/WH/008_AM/ 20111020_1297_8ch_1mm3/output/OutPutFiles_1/BrainMaps/WF_map.nii.gz');
BMs1=readFileNifti('/home/data/Stanford_Data/WH/008_AM/ 20111020_1297_8ch_1mm3/output/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz');
T1S1=readFileNifti('/home/data/Stanford_Data/WH/008_AM/ 20111020_1297_8ch_1mm3/output/OutPutFiles_1/T1w/T1w1.nii.gz');
T1S1.data=10.*T1S1.data;
% csf32=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1208_32ch_1mm/output_WLM0/SPGR_1/Align_0.9375_0.9375_1/csf_seg_T1_large.nii.gz');

WFS2=readFileNifti('/home/data/Stanford_Data/WH/008_AM/20111020_1208_32ch_1mm/output/OutPutFiles_1/BrainMaps/WF_map.nii.gz');
BMS2=readFileNifti('/home/data/Stanford_Data/WH/008_AM/20111020_1208_32ch_1mm/output/SPGR_1/Align_0.9375_0.9375_1/brainMask.nii.gz');
T1S2=readFileNifti('/home/data/Stanford_Data/WH/008_AM/20111020_1208_32ch_1mm/output/OutPutFiles_1/T1w/T1w1.nii.gz');
T1S2.data=10.*T1S2.data;
% csf8=readFileNifti('/home/shai.berman/Documents/Code/AFQ_test/008_AM/20111020_1297_8ch_1mm3/output_WLM0/SPGR_1/Align_0.9375_0.9375_1/csf_seg_T1_large.nii.gz');


BM = logical(BMS2.data) & logical(BMs1.data);% & (~logical(csf32.data)) & (~logical(csf8.data));

% This region is cut off at the edge of scan box
BM(:,:,1:52)   = 0;
BM(:,:,144:end)= 0;

% should we do it with T1. We want T1 that is greater then 0.5 and smaller
% then 2.5.
BM=BM & T1S1.data>0.5 & T1S1.data<2.5 & T1S2.data>0.5 & T1S2.data<2.5;

% BM=BM & WF8ch.data<0.95 & WF32ch.data<0.95;

% The coefficient of determination (R^2)

CV = (calccod(WFS2.data(BM),WFS1.data(BM))/100);


%% Plot a  2D histogram of the error

% We plot the MTV = 1 - WF the macromolecular tisuue volume fraction.
nBins = 155;
% [n,x,y] = mrQ_hist2d(WF8ch.data(BM),WF32ch.data(BM),nBins);
[n,x,y] = mrQ_hist2d(WFS2.data(BM),WFS1.data(BM),nBins);

maxN = ceil(max(n(:))/10)*10;


mrvNewGraphWin;
image(x(1,:),y(:,1),uint8(n./maxN.*255));
colormap(flipud(gray(256)));
identityLine(gca);
axis square xy; axis image;
xlim([0.6 1]); ylim([0.6 1])
ylabel('PD  8ch' ,'FontSize',16);
xlabel('PD  32ch' ,'FontSize',16);
set(gca,'FontSize',16)
title( [' R^2= ' num2str(CV) '  WM0'])
grid on

showMontage(WFS2.data(:,:,90)); caxis([0 1]);axis off
showMontage(WF32chWarp.data(:,:,90));caxis([0 1]);axis off



