function mrQ=mrQ_correctB1(mrQ,B1fileName, T1app,T1_GS,outDir)



%% I. Check INPUTS and set defaults

 if notDefined('outDir');
    outDir =mrQ.Ants_Info.SEIR_SPGR_Curent_AlignDirs{2};
end

if notDefined('B1fileName')
    B1fileName=mrQ.B1.unCorrected;
end

if notDefined('T1app')
   T1app=mrQ.T1_B1_LFit_unCorrected;           
end

if notDefined('BMapp')
   BMapp=mrQ.LinFit.BrainMask;
end
if notDefined('T1_GS')
    T1_GS=mrQ.SEIR_epi_T1file;
end

if notDefined('BM_GS')
    BM_GS=mrQ.SEIR_epi_Maskfile;
end

%% find the calibration value of B1
    % load the gold standard SEIR
    seirT1=readFileNifti(T1_GS);
    seirBM=readFileNifti(BM_GS);
    seirBM=logical(seirBM.data) & seirT1.data<1000 & seirT1.data>700;
    
    % load the "apparent" T1 - SPGR
    spgrT1=readFileNifti(T1app);
    spgrBM=readFileNifti(BMapp);
    spgrBM=logical(spgrBM.data) & spgrT1.data<1 & spgrT1.data>0.7;
    
    % find the maxima around the WM values
    [seir_Dens, seir_Vals] = ksdensity(seirT1.data(seirBM));
    mxT1_seir = ( seir_Vals(seir_Dens==max(seir_Dens)) ) / 1000;
    
    [spgr_Dens, spgr_Vals] = ksdensity(spgrT1.data(spgrBM));
    mxT1_spgr = spgr_Vals(spgr_Dens==max(spgr_Dens));
    
    % find the ratio between the two peaks
    sc = mxT1_spgr/mxT1_seir;
    
    
    %% load and correct uncorrected B1
    
    B1=readFileNifti(B1fileName);
    B1.data=(sqrt(sc))*B1.data;
    
    %% save the new B1
   B1filePath=fullfile(outDir,'B1_Map.nii.gz');
   mrQ.B1FileName=B1filePath;
   dtiWriteNiftiWrapper(B1.data,B1.qto_xyz, B1filePath);
   save(mrQ.name,'mrQ')
   fprintf(['Done fitting B1 map. B1 file is saved: '   mrQ.B1FileName        '  \n']);
