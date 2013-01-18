function mrQ_FIT_MT(mrQ,mtInd,t1,m0,s0,tt1,B1)

load(fullfile(mrQ.MT.MT_initDir,'dat_aligned.mat'));

if notDefined('mtInd')
    mtInd=(ones(1,length(s)));
end

if notDefined('t1')
t1file=fullfile(mrQ.spgr_initDir,'maps/T1_map_lsq.nii.gz');%
t1=readFileNifti(t1file);t1=t1.data;
end


if notDefined('B1')
B1file=fullfile(mrQ.spgr_initDir,'B1_Map.nii.gz');%
B1=readFileNifti(B1file);B1=B1.data;
end



fitMethod = 'ls';
mtInd=logical(mtInd);
mtFlip = s(mtInd(1)).flipAngle;
mtTr = s(mtInd(1)).TR;
mtData = cat(4,s(mtInd).imData);

brainMaskfile=fullfile(mrQ.spgr_initDir,'brainMask');%
if exist(file,brainMaskfile)
    brainMask=readFileNifti(brainMaskfile);brainMask=brainMask.data(:);
else
    brainMask=find(mtData(:,:,:,1));
end

delta=[s(mtInd).mtOffset];


if notDefined('s0')
    if notDefined('m0')
        m0file=fullfile(mrQ.spgr_initDir,'M0_LFit_HM.nii.gz');%
        m0=readFileNifti(m0file);m0=m0.data;
    end
    if notDefined('tt1')
        tt1file=fullfile(mrQ.spgr_initDir,'T1_LFit_HM.nii.gz');%
        tt1=readFileNifti(tt1file);tt1=tt1.data;
    end
    fa = mtFlip*pi/180;
    
    s0 = m0 .* sin(fa) .* (1-exp(-mtTr./tt1)) ./ (1-cos(fa) .* exp(-mtTr./tt1));
end


[f,k,s0,gof] = mrQ_relaxFitMtDist(mtData, delta, s0, pd, mtTr, mtFlip, brainMask, fitMethod, B1);
