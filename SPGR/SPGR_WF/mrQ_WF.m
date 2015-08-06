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


if notDefined('PDfile')
    load(mrQ.opt_logname);
    PDfile=opt.PDfile;
end

% 
% M0=readFileNifti(opt.M0file_Org);outMm=M0.pixdim(1:3);     M0=M0.data; MSZ=size(M0);
%     %M0=sqrt(sum(M0.^2,4));
%     M0=(sum(M0,4));


%%
%1 find area of high T1=4.3 +-0.1 around the vetrical.
%2. loop over raw (align. and compbined) data and find the flip angle(s)
%with the SNR. it is hard to measure noise. so we will use the signal in WM as a refernce becouse this is a place with good signal across flip angles) 
%3.exclude M0 outliers in the ROI
%4. calculate  PD within the selected ROI and flip angle(s) (using the SPGR eqation, T1, and gain [note gain for multi coil is define for the sum of multi coils, and for the combine as arrived from the scanner]) 
%5 scale this roi to have pd of 1. --> this is our gloobl scale
%6. apllay to PD images to make WFfile --> done
