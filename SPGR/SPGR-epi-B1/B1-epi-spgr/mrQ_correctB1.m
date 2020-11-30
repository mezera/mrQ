function mrQ=mrQ_correctB1(mrQ,B1fileName, T1app,T1_GS,outDir)



%% I. Check INPUTS and set defaults

if notDefined('outDir')
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

if ~isfield(mrQ,'ShiftB1')
    mrQ.ShiftB1=0;
end
%% find the calibration value of B1

% The calibration depand on the field strength (S.F):
if ~isfield(mrQ,'fieldstrength')
    error('Uknown field strength')
end
    if mrQ.fieldstrength==3
        lb=700;
        ub=1000;
    elseif mrQ.fieldstrength==1.5
        lb=650;
        ub=900;
    else
        error('Uknown field strength')
    end
% load the gold standard SEIR
    seirT1=readFileNifti(T1_GS);
    seirBM=readFileNifti(BM_GS);
    seirBM=logical(seirBM.data) & seirT1.data<ub & seirT1.data>lb;
    %      seirBM=logical(seirBM.data) & seirT1.data<1500 & seirT1.data>700;
    
    % load the "apparent" T1 - SPGR
    spgrT1=readFileNifti(T1app);
    spgrBM=readFileNifti(BMapp);
    spgrBM=logical(spgrBM.data) & spgrT1.data<ub/1000 & spgrT1.data>lb/1000;
  %  spgrBM=logical(spgrBM.data) & spgrT1.data<1.5 & spgrT1.data>0.7;
    
    % find the maxima around the WM values
    [seir_Dens, seir_Vals] = ksdensity(seirT1.data(seirBM));
    mxT1_seir = ( seir_Vals(seir_Dens==max(seir_Dens)) ) / 1000;
    
    [spgr_Dens, spgr_Vals] = ksdensity(spgrT1.data(spgrBM));
    mxT1_spgr = spgr_Vals(spgr_Dens==max(spgr_Dens));
    
    % find the ratio between the two peaks
    sc = mxT1_spgr/mxT1_seir;
    
    mrQ.B1shiftpeaks=[mxT1_spgr mxT1_seir];
    save(mrQ.name,'mrQ')

    if (sc<0.9 || sc>1.1)  && mrQ.ShiftB1==0
        errormassage=['There is a very big shift between the white matter peak in T1-spgr and the white matter peak in T1-seir.\n' ...
            'there might be a problem with the way we find these peaks.\n' ...
            'It could be related to the seir brain mask, or to the data itself.\n' ...
            'Please check that both T1 in SPGR space and T1 in SEIR space have a clear white matter peak and that the seir brain mask is good. \n' ...
            'The T1-spgr path: %s \n' ...
            'The T1-seir path: %s \n'...
            'The brain mask-seir path: %s \n'...
            'Please check the field mrQ.B1shiftpeaks=[T1_spgr_peak T1_seir_peak] ' ...
            'and manually test that the peaks were found correctly. \n ' ...
            'An additional B1 adjustment is implemented in the code to shift the spgr values toward the seir values, \n' ...
            'However, in this case this adjustment is not recommended without further validation of the data '];
        error( 'u:stuffed:it' , errormassage,mrQ.T1_B1_LFit_unCorrected,mrQ.SEIR_epi_T1file,mrQ.SEIR_epi_Maskfile)
    end
    if (sc<0.95 || sc>1.05) && mrQ.ShiftB1==0
        errormassage=['The white matter peak in T1-spgr is slightly shifted from the white matter peak in T1-seir. \n ' ...
            'This shift might require additional adjustment to the B1 map. \n ' ...
            'Please check the field mrQ.B1shiftpeaks=[T1_spgr_peak T1_seir_peak] ' ...
            'and manually test that the peaks were found correctly. \n ' ...
            'Check the T1-spgr: %s \n' ...
            'The T1-seir: %s \n'...
            'and the brain mask-seir path: %s \n'...
            'To implement the additional B1 adjustment please ' ...
            'manually change mrQ.ShiftB1 to be =1, save mrQ, and rerun.\n' ...
            'To proceed without further adjustment of B1, manually change mrQ.ShiftB1 to be =2, save mrQ, and rerun.' ];

        error( 'u:stuffed:it' , errormassage,mrQ.T1_B1_LFit_unCorrected,mrQ.SEIR_epi_T1file,mrQ.SEIR_epi_Maskfile)
    end
    %% load and correct uncorrected B1


    B1=readFileNifti(B1fileName);
    if mrQ.ShiftB1==1
        B1.data=(sqrt(sc))*B1.data;
        fprintf(['Additional B1 adjustment (peak shift) was done \n']);
    end
    %% save the new B1
   B1filePath=fullfile(outDir,'B1_Map.nii.gz');
   mrQ.B1FileName=B1filePath;
   mrQ.B1shiftpeaks=[mxT1_spgr mxT1_seir];
   dtiWriteNiftiWrapper(B1.data,B1.qto_xyz, B1filePath);
   save(mrQ.name,'mrQ')
   fprintf(['Done fitting B1 map. B1 file is saved: '   mrQ.B1FileName        '  \n']);
