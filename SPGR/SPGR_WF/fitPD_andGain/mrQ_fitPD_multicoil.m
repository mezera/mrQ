function [opt]=mrQ_fitPD_multicoil(outDir,SunGrid,M0cfile,degrees,sqrtF,subName,proclus,clobber)
%
% mrQ_fitPD_multicoil(outDir,SunGrid,M0cfile,degrees,sqrtF)
% # Create the PD from the boxes
%
% INPUTS:
%   outDir      - The output directory - also reading file from there
%
%   SunGrid     - Flag to use the SGE for computations
%   M0cfile     - The combined/aligned M0 data
%   degrees     - Polynomial degrees for the coil estimation
%   sqrtF       - use the sqrt of the signal (defult is no zero)
%   subName     - the subject name for the SGE run
%  clobber:     - Overwrite existing data and reprocess. [default = false]
% OUTPUTS:
% opt           - a  structure that save the fit parameter is saved in the outdir and in
%               the tmp directory named fitLog.mat
%               calling to FitM0_sanGrid_v2 that save the output fitted files in tmp directorry
%              this will be used lster by mrQfitPD_multiCoils_M0 to make the PD map
%
%
% SEE ALSO:
%   mrQmulticoilM0.m
%
%
% ASSUMPTIONS:
%   This code assumes that your SGE output directory is '~/sgeoutput'
%
%
% (C) Stanford University, VISTA
%
%


%% CHECK INPUTS AND SET DEFAULTS

if (notDefined('outDir') || ~exist(outDir,'dir'))
    outDir = uigetDir(pwd,'Select outDir');
end

if(~exist('degrees','var'))
    disp('Using the defult polynomials: Degrees = 2 for coil estimation');
    degrees = 3;
end

% if the M0cfile is not an input  we load the file that was made by
% mrQ_multicoilM0.m (defult)
if(~exist('M0cfile','var') || isempty(M0cfile))
    M0cfile{1} = fullfile(outDir,'AligncombineCoilsM0.nii.gz');
    
    if ~exist(M0cfile{1},'file')
        disp(' can not find the multi coils M0 file')
        error
    end
end
%using sungrid  defult no
if(~exist('SunGrid','var'))
    disp('SGE will not be used.');
    SunGrid = false;
end

% sqrt the M0 file we can do it and undo it in the end (defult is not to do
% it)
if(~exist('sqrtF','var')  || isempty(sqrtF))
    sqrtF=0;
end

% Clobber flag. Overwrite existing fit if it already exists and redo the PD
% Gain fits
if notDefined('clobber')
    clobber = false;
end

%%

% we set an output strcture opt that will send to the grid with all the
% relevant information for the Gain fit
opt{1}.degrees = degrees;
opt{1}.outDir  = outDir;
opt{1}.sqrtF=sqrtF;

% Load the brain mask file from outDir
BMfile = fullfile(outDir,'brainMask.nii.gz');

if (exist(BMfile,'file'))
    disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    % xform   = brainMask.qto_xyz;
    mmPerVox  = brainMask.pixdim;
    brainMask = logical(brainMask.data);
else
    error('Cannot find the file: %s', BMfile);
end


% Load the CSF file from mrQ_CSF.m
cdffileFS = fullfile(outDir,'csf_seg_T1.nii.gz');
if (exist(cdffileFS,'file'))
    M0f = readFileNifti(cdffileFS);
    opt{1}.M0f = double(M0f.data);
else
    error('Cannot find CSF file: %s',cdffileFS);
end

% Get the subject prefix for SGE job naming
if notDefined('subName')
    % This is a job name we get from the for SGE
    [~, subName] = fileparts(fileparts(fileparts(fileparts(fileparts(outDir)))));
    disp([' Subject name for lsq fit is ' subName]);
end





%% Define the mask that will be used to find the boxes to fit

Avmap = zeros(size(brainMask));


% Create morphological structuring element (STREL)
NHOOD         = [0 1 0; 1 1 1; 0 1 0];
NHOOD1(:,:,2) = NHOOD;
NHOOD1(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
NHOOD1(:,:,3) = NHOOD1(:,:,1);
SE            = strel('arbitrary', NHOOD1);

% Dilate the structuring element
% bb = imdilate(b,SE);

% Try to fill the brain with boxes of roughly 20mm^3 and with an
% overlap of 2;
sz   = (size(brainMask));
boxS = round(24./mmPerVox);
even = find(mod(boxS,2)==0);

boxS(even)  = boxS(even)+1;
opt{1}.boxS = boxS;

% Determine the percentage of overlap (pol) (0.1, 0.5, 0.7)
pol     = 0.5;
overlap = round(boxS.*pol);

% The center of the boxs we will fit
[opt{1}.X,opt{1}.Y,opt{1}.Z] = meshgrid(round(boxS(1)./2):boxS(1)-overlap(1):sz(1)    ,    round(boxS(2)./2):boxS(2)-overlap(2):sz(2)    , round(boxS(3)./2):boxS(3)-overlap(3):sz(3));

controlmask = zeros(size(opt{1}.X));
donemask    = controlmask;

opt{1}.HboxS = (boxS-1)/2;

% Two masks that control what box done and which need to e done
% box = CSF(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3));



%% Intiate parameters for the fit
%options =  optimset('LevenbergMarquardt','on','Tolx',1e-12,'TolF',1e-12,'MaxIter',100,'MaxFunEvals',inf,'Display', 'iter');%
a=version('-date');
if str2num(a(end-3:end))==2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'off','Tolx',1e-12);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'off','Tolx',1e-12);%'TolF',1e-12
    
end

opt{1}.numIn = 8; %number of coil we will use to fit (best 8);
% [opt{1}.Poly,opt{1}.str] = constructpolynomialmatrix3d(boxS,find(ones(boxS)),degrees);
% opt{1}.lb = ones(opt{1}.numIn,size(opt{1}.Poly,2)).*-inf;
% opt{1}.ub = ones(opt{1}.numIn,size(opt{1}.Poly,2)).*inf;




%% Loop over the box to fit and check there is data there

ii = 1;
opt{1}.donemask = donemask;

for i=1:prod(size(opt{1}.X))
    
    [fb(1) fb(2) fb(3)] = ind2sub(size(opt{1}.X),i);
    [empty,Avmap] = mrQM0fit_isBox(opt{1},brainMask,fb,Avmap);
    
    if empty == 0
        opt{1}.wh(ii) = i;
        ii = ii+1;
        
    elseif empty == 1
        donemask(fb(1),fb(2),fb(3)) = -1e3;
        
    elseif empty == -1
        donemask(fb(1),fb(2),fb(3)) = -2e3;
        opt{1}.wh(ii) = i;
        ii = ii+1;
    end
    
end

opt{1}.donemask = donemask;


%% Loop over the combined M0 files and perform the fits for each box

% tipicaly we will combine the raw Mo images first in mrQ_multicoilM0 so this will be 1 file

for j=1:length(M0cfile)
    
    % Read the M0c image
    %  M = readFileNifti(M0cfile{j});
    %     opt{1}.coils = size(M.data,4);
    
    sgename    = [subName '_MultiCoilM0' num2str(j)];
    dirname    = [outDir '/tmpSGM0' num2str(j)];
    dirDatname = [outDir '/tmpSGM0dat' num2str(j)];
    jumpindex  = 10; %number of boxs fro each SGR run
    
    opt{1}.dat = M0cfile{j};
    opt{1}.dirDatname = dirDatname;
    opt{1}.name = [dirname '/M0boxfit_iter'] ;
    opt{1}.date = date;
    opt{1}.jumpindex = jumpindex;
    opt{1}.brainMask = brainMask;
    opt{1}.SGE=sgename;
    % Save out a logfile with all the options used during processing
    logname = [outDir '/fitLog.mat'];
    
    save(logname,'opt');
    
    
    if clobber && (exist(dirname,'dir'))
        % in the case we start over and there are  old fits, so we will
        % deleat them
        eval(['! rm -r ' dirname]);
    end
    clear  brainMask donemask options controlmask M0f
    %%   Perform the gain fits
    % Perform the fits for each box using the Sun Grid Engine
    if SunGrid==1;
        
        % Check to see if there is an existing SGE job that can be
        % restarted. If not start the job, if yes prompt the user.
        if (~exist(dirname,'dir')),
            mkdir(dirname);
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            if proclus==1
                sgerun2('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',sgename ,1,1:ceil(length(opt{1}.wh)/jumpindex),[],[],8000); 
            else
                sgerun('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt{1}.wh)/jumpindex),[],[],8000);
            end
        else
            % Prompt the user
            inputstr = sprintf('An existing SGE run was found. \n Would you like to try and finish the exsist SGE run?');
            an1 = questdlg( inputstr,'mrQ_fitPD_multiCoils','Yes','No','Yes' );
            if strcmpi(an1,'yes'), an1 = 1; end
            if strcmpi(an1,'no'),  an1 = 0; end
            
            % User opted to try to finish the started SGE run
            if an1==1
                reval = [];
                list  = ls(dirname);
                ch    = 1:jumpindex:length(opt{1}.wh);
                k     = 0;
                
                for ii=1:length(ch),
                    ex=['_' num2str(ch(ii)) '_'];
                    if length(regexp(list, ex))==0,
                        k=k+1;
                        reval(k)=(ii);
                    end
                end
                
                if length(find(reval)) > 0
                    eval(['!rm -f ~/sgeoutput/*' sgename '*'])
                    if proclus==1
                        for kk=1:length(reval)
                        sgerun2('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',[sgename num2str(kk)],1,reval(kk),[],[],3000);
                        %sgerun2('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',[sgename 'redo'],1,reval,[],[],8000);

                        end
                    else
                        sgerun('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',sgename,1,reval,[],[],8000);
                    end
                end
                
                % User opted to restart the existing SGE run
            elseif an1==0,
                t = pwd;
                cd (outDir)
                eval(['!rm -rf ' dirname]);
                cd (t);
                eval(['!rm -f ~/sgeoutput/*' sgename '*'])
                mkdir(dirname);
                if proclus==1
                    sgerun2('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt{1}.wh)/jumpindex),[],[],8000);
                else
                    sgerun('FitM0_sanGrid_v2(opt,jumpindex,jobindex);',sgename,1,1:ceil(length(opt{1}.wh)/jumpindex),[],[],8000);
                    
                end
            else
                error('User cancelled');
            end
        end
    else 
         fprintf('\n fit the PD map localy, may be slow. SunGrid use can be much faster             \n');

        if (~exist(dirname,'dir')),
            mkdir(dirname);
            jobindex=1:ceil(length(opt{1}.wh)/jumpindex);
        else
              jobindex = [];
                list  = ls(dirname);
                ch    = 1:jumpindex:length(opt{1}.wh);
                k     = 0;
                
                for ii=1:length(ch),
                    ex=['_' num2str(ch(ii)) '_'];
                    if length(regexp(list, ex))==0,
                        k=k+1;
                        jobindex(k)=(ii);
                    end
                end
                
            
        end
        
         if ~isempty(jobindex)
        for i=jobindex
            FitM0_sanGrid_v2(opt,jumpindex,i);
        end
        end
       
    end
    
end

