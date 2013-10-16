function mrQ_antsCallBypassScript(outDir,SEIRepi_Dir,AnalysisInfo,mrQ.AligndSPGR)



%%  make the files  names
%t1fileHM = fullfile(outDir,['T1_LFit_HM.nii.gz']);
 %SET1file=[SEIRepi_Dir '/fitT1_GS/T1FitNLSPR_SEIR_Dat_T1.nii.gz'];
AnalysisInfo.T1_spgr_epi_RB=fullfile(outDir,'Warp_T1_SPGRT2EPI_RB.nii.gz');


for d=1:length(flipAngles) %loof over
     AnalysisInfo.Raw_spgr_epi_RB{d}=[AligndSPGR{d} 'WarpRB_SPGRT2EPI.nii.gz'];
      %% save and document 

 end 

morefiles2_HR=mrQ.AligndSPGR;
 
 InsavefileN{1}=AnalysisInfo.T1_spgr_epi_RB;

for d=1:length(flipAngles) %loof over
    InsavefileN{d+1}=AnalysisInfo.Raw_spgr_epi_RB{d};

 end 
AnalysisInfo.ManualWarp=1;
%%
%Align manually
mrQ_compare2Mpas(SET1file,t1fileHM,outDir,morefiles2_HR,InsavefileN)
%mrQ_knkREgisterIm(SET1file,t1fileHM,KNK_ManParam,outDir,morefiles2_HR,InsavefileN);


%%

%make a structuure to work with for B1 fit
t1seir=readFileNifti(SET1file);
Res{1}.im=t1seir.data;
Res{1}.name='t1SEIRepi';
clear t1seir
t1spgr=readFileNifti(AnalysisInfo.T1_spgr_epi_RB);
Res{2}.im=t1spgr.data;
Res{2}.name='t1SPGR_in_epi_space';

for i=1:length(flipAngles)
  im=readFileNifti(AnalysisInfo.Raw_spgr_epi_RB{i}); 
Res{i+2}.im=im.data;
Res{i+2}.name=['align_rawFA' num2str(flipAngles(i))] ;

end;
%save the strucute

  save(AlignFile,'Res');

   %document  and done
AnalysisInfo.spgr2epi_Align_date=date;

infofile=fullfile(outDir,'AnalysisInfo.mat');
save(infofile,'AnalysisInfo');
  
    

