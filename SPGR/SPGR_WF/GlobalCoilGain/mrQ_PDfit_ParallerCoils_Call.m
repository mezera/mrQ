function [opt]=mrQ_PDfit_ParallerCoils_Call(outDir,SunGrid,M0cfile,degrees,subName,proclass,prctileClip,Ncoils,Nfits,clobber)
%
% [opt]=mrQ_PDfit_ParallerCoils(outDir,SunGrid,M0cfile,degrees,subName,proclass,clobber)
% # Create the PD from the boxes
%
% INPUTS:
%   outDir      - The output directory - also reading file from there
%
%   SunGrid     - Flag to use the SGE for computations
%   M0cfile     - The combined/aligned M0 data
%   degrees     - Polynomial degrees for the coil estimation
%   subName     - the subject name for the SGE run
%   prctileClip - the part of the data we clip becouse it is changing too
%                  fast in space and hard to be fit by poliynomyials
%   Ncoils      - number of coil use for  each fit  (defult 10)
%   Nvalus      - number of voxel use for  each fit  (defult 10000)
%   Nfits       - number of fits we will run (and avrage after)
%  clobber:     - Overwrite existing data and reprocess. [default = false]
% OUTPUTS:
% opt           - a  structure that save the fit parameter is saved in the outdir and in
%               the tmp directory named fitLog.mat
%               calling to FitM0_sanGrid_v2 that save the output fitted files in tmp directorry
%              this will be used lster by mrQfitPD_multiCoils_M0 to make the PD map
%
%
% SEE ALSO:
%   
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

if notDefined('degrees')
    disp('Using the defult polynomials: Degrees = 3 for coil estimation');
    degrees = 3; 
end

%using sungrid  defult no
if(~exist('SunGrid','var'))
    disp('SGE will not be used.');
    SunGrid = false;
end


if notDefined('Ncoils')
    Ncoils = 10;
end
if notDefined('Nvalus')
    Nvalus = 10000;
end
if notDefined('Nfits')
    Nfits = 200;
end

% Clobber flag. Overwrite existing fit if it already exists and redo the PD
% Gain fits
if notDefined('clobber')
    clobber = false;
end

% we set an output strcture opt that will send to the grid with all the
% relevant information for the Gain fit
opt{1}.degrees = degrees;
opt{1}.outDir  = outDir;
opt{1}.Ncoils  = Ncoils;
opt{1}.Nvalus  = Nvalus;
opt{1}.Nfits  = Nfits;


% Get the subject prefix for SGE job naming
if notDefined('subName')
    % This is a job name we get from the for SGE
    [~, subName] = fileparts(fileparts(fileparts(fileparts(fileparts(outDir)))));
    disp([' Subject name for lsq fit is ' subName]);
end


sgename    = [subName '_MultiCoilM0'];
dirname    = [outDir '/tmpSGM0_P3' ];
opt{1}.name = [dirname '/M0boxfit_iter'] ;
opt{1}.date = date;

opt{1}.SGE=sgename;
    % Save out a logfile with all the options used during processing
    
    %% select the coil we will work on
    BMfile = fullfile(outDir,'HeadMask.nii.gz');
    if (exist(BMfile,'file'))
       
    else
        error('Cannot find the file: %s', BMfile);
    end
   opt{1}.BMfile=BMfile;
    % if the M0cfile is not an input  we load the file that was made by
    % mrQ_multicoilM0.m (defult)
    if(~exist('M0cfile','var') || isempty(M0cfile))
        M0cfile = fullfile(outDir,'AligncombineCoilsM0.nii.gz');
    end
    if ~exist(M0cfile,'file')
        disp(' can not find the multi coils M0 file')
        error
   
    end

 opt{1}.M0cfile = M0cfile;

% we need to find the coils data that is usful (couple of coils that can be
% mached). and define the area of the coil that can be fitted by
% polynoyials ((above the noise floor and not to varing to fast in space).
if notDefined('prctileClip')
    prctileClip=[];
end
 opt{1}.prctileClip =prctileClip;

%%
TMfile=fullfile(outDir,'FS_tissue.nii.gz');
 if ~exist(TMfile,'file')
        disp(' can not find the multi coils segmetation file')
        error
        
    end
opt{1}.TM=TMfile;



%% send it to be fitted on the grid

logname = [outDir '/PDGridCall_Log.mat'];
    
    save(logname,'opt');
    if clobber && (exist(dirname,'dir'))
        % in the case we start over and there are  old fits, so we will
        % deleat them
        eval(['! rm -r ' dirname]);
    end
    %%   Perform the gain fits
    % Perform the fits for each box using the Sun Grid Engine
    if SunGrid==1;
        
        % Check to see if there is an existing SGE job that can be
        % restarted. If not start the job, if yes prompt the user.
        if (~exist(dirname,'dir')),
            mkdir(dirname);
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            if proclass==1
                sgerun2('mrQ_PDrandfit_ParallerCoils_Gridcall(logname,jobindex);',sgename,1,1:Nfits,[],[],15000);

            else 
                sgerun('mrQ_PDrandfit_ParallerCoils_Gridcall(logname,jobindex);',sgename,1,1:Nfits,[],[],15000);
                
            end
        else
            % Prompt the user
            inputstr = sprintf('An existing SGE run was found. \n Would you like to try and finish the exsist SGE run?');
            an1 = questdlg( inputstr,'mrQ_PDfit_ParallerCoils_Call','Yes','No','Yes' );
            if strcmpi(an1,'yes'), an1 = 1; end
            if strcmpi(an1,'no'),  an1 = 0; end
            
            % User opted to try to finish the started SGE run
            if an1==1
                reval = [];
                list  = ls(dirname);
                ch    = 1:Nfits;
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
                    if proclass==1
                        for kk=1:length(reval)
                            sgerun2('mrQ_PDrandfit_ParallerCoils_Gridcall(logname,jobindex);',[sgename num2str(kk)],1,reval(kk),[],[],15000);

                        end
                    else
                            sgerun('mrQ_PDrandfit_ParallerCoils_Gridcall(logname,jobindex);',sgename ,1,reval,[],[],15000);
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
                if proclass==1
                sgerun2('mrQ_PDrandfit_ParallerCoils_Gridcall(logname,jobindex);',sgename,1,1:Nfits,[],[],15000);
                else
                sgerun('mrQ_PDrandfit_ParallerCoils_Gridcall(logname,jobindex);',sgename,1,1:Nfits,[],[],15000);
                    
                end
            else
                error('User cancelled');
            end
        end
        
    end
    



