%Fig_PD_G_and_M0_examples
IM2=readFileNifti('/biac2/wandell2/data/WH/001_RM/Qmr/20120106_1648/SPGR_1/Align_0.9375_0.9375_1/maps/WF_map.nii.gz');
figure;imagesc(IM2.data(:,:,90));colormap gray
mrUtilResizeFigure(gcf, 900, 900);
mrUtilPrintFigure('PD_example_slice.eps');

clear IM2
IMM0=readFileNifti('/biac2/wandell2/data/WH/001_RM/Qmr/20120106_1648/SPGR_1/Align_0.9375_0.9375_1/AligncombineCoilsM0.nii.gz');
coils=[4 6 18 28];
for ii=coils
    figure;imagesc(IMM0.data(:,:,90,ii));colormap gray;axis off
    
mrUtilResizeFigure(gcf, 900, 900);
keyboard
mrUtilPrintFigure(['M0_example_slice_coil_' num2str(ii) '.eps']);
end
clear IMM0
IMG=readFileNifti('/biac2/wandell2/data/WH/001_RM/Qmr/20120106_1648/SPGR_1/Align_0.9375_0.9375_1/Gains.nii.gz');


for ii=coils
    figure;imagesc(IMG.data(:,:,90,ii));colormap gray;axis off
    
mrUtilResizeFigure(gcf, 900, 900);

mrUtilPrintFigure(['CoilG_example_slice_coil_' num2str(ii) '.eps']);
end