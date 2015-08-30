function mrQ=mrQ_arrangeOutPutDir(mrQ)
%mrQ=mrQ_arrangeOutPutDir(mrQ)
%
% This function rearranges the data from the mrQ structure, as it organizes
% the previous created NIfTI files into a convenient directory called
% 'OutPutFiles'. This contains three subdirectories:
%            1. BiasMap
%                 (files: Gains, B1Map)
%            2. T1w
%                 (files: T1w)
%            3. BrainMaps
%                 (files: SIR, T1, TV, VIP)
%                 SEIR
%                      (files: T1)
%
% This function accepts the mrQ structure as its INPUT, and returns the
% updated mrQ structure as its OUTPUT.
%
%

%% I. Create an Output directory in the mrQ structure

fprintf('\n Organize the output directory                \n');

upDir=fileparts(mrQ.name);
OutPutNiiDirName=fullfile(upDir,'OutPutFiles');

%if (exist(OutPutNiiDir,'dir'))
ex=false; num=1;

while ex==false
    OutPutNiiDir=[OutPutNiiDirName '_' num2str(num)] ;
    
    if (~exist(OutPutNiiDir,'dir'))
        mkdir(OutPutNiiDir);
        mrQ.OutPutNiiDir=OutPutNiiDir;
        ex=true;
    else
        num=num+1;
    end
    
end

%% II. Organize the brain maps in a directory named 'BrainMaps'
fprintf('\n Organized the maps in a maps directory                \n');
mapDir=fullfile(mrQ.OutPutNiiDir,'BrainMaps');

% If we redo it and the maps directory already exists,
%   we will save the old one before we make a new one.

if (exist(mapDir,'dir'))
    ex=0; num=0;
    while ex==0
        num=num+1;
        mapDirOld=fullfile(mrQ.OutPutNiiDir,['mapsOld_' num2str(num)] );
        if (~exist(mapDirOld,'dir'))

            eval(['! mv ' mapDir  ' ' mapDirOld]);
            ex=1;
        end
    end
end

mkdir(mapDir);
[T1file, M0file,BMfile]=mrQ_get_T1M0_files(mrQ,1,1,1);

  [d dd ddd]=fileparts(T1file);
 cmd =(['! mv ' T1file ' ' mapDir '/.']);
     eval(cmd);
  
    cmd =(['! ln -s  '  mapDir '/' dd ddd ' ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

% % If T1_map_lsq.nii.gz exists in the spgr_intDir, then move it to the map
% % directory. Otherwise, assume that it ended up in mapDirOld and retain
% % this T1 map. This would occur if the PD fits were re-run but the T1
% % fits were not.
% %
% if exist(fullfile(mrQ.spgr_initDir, 'T1_map_lsq.nii.gz'),'file')
%     cmd =(['! mv ' mrQ.spgr_initDir '/T1_map_lsq.nii.gz ' mapDir '/.']);
%     eval(cmd);
%     
%     cmd =(['! ln -s  '  mapDir '/T1_map_lsq.nii.gz ' mrQ.spgr_initDir '/.']) ;
% eval(cmd);
% 
% else
%     % If a new T1 map was not generated, then copy it from the old
%     % directory. We copy rather than move to leave the old directory intact
%     %
%     cmd =(['! cp ' mapDirOld '/T1_map_lsq.nii.gz ' mapDir '/.']);
%     eval(cmd);
% end
%% need to add a check if the maps exist in the mrQ.spgr_initDir or in an outputfiles dir

cmd =(['! mv ' mrQ.spgr_initDir '/WF_map.nii.gz ' mapDir '/.']) ; eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/VIP_map.nii.gz ' mapDir '/.']) ;eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/TV_map.nii.gz ' mapDir '/.']) ;eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/SIR_map.nii.gz ' mapDir '/.']) ;eval(cmd);

mrQ.mapsDir=mapDir;
%Which T1 need a fix
mrQ.maps.T1path=fullfile(mapDir,[dd, ddd]);
%% 

mrQ.maps.WFpath=fullfile(mapDir,'WF_map.nii.gz');
mrQ.maps.TVpath=fullfile(mapDir,'TV_map.nii.gz');
mrQ.maps.SIRpath=fullfile(mapDir,'SIR_map.nii.gz');
mrQ.maps.VIPpath=fullfile(mapDir,'VIP_map.nii.gz');

%% II-a. Copy SEIR images to a subdirectory within 'Brain Maps'

SEIRT1Dir=fullfile(mapDir,'SEIR');
mkdir(SEIRT1Dir)
T1SEIRfile=mrQ.SEIR_epi_T1file;
% cmd =(['! cp ' mrQ.T1file ' '  SEIRT1Dir '/.']) ;
cmd =(['! cp ' T1SEIRfile ' '  SEIRT1Dir '/.']) ;
eval(cmd);

%% III. Organize the bias maps in a directory named 'BiasMaps'

BiasDir=fullfile(mrQ.OutPutNiiDir,'BiasMap');

% If we redo it and maps directory already exists,
%   we will save the old one before we make a new one.
if (exist(BiasDir,'dir'))
    ex=0; num=0;
    while ex==0
        num=num+1;
        BiasDirOld=fullfile(mrQ.OutPutNiiDir,['BiasDirOld_' num2str(num)] );
        if (~exist(BiasDirOld,'dir'))

            eval(['! mv ' BiasDir  ' ' BiasDirOld]);
            ex=1;
        end
    end
end

mkdir(BiasDir);

%
cmd =(['! mv ' mrQ.spgr_initDir '/B1_Map.nii.gz ' BiasDir '/.']) ;
eval(cmd);

cmd =(['! ln -s  '  BiasDir '/B1_Map.nii.gz ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/Gains.nii.gz ' BiasDir '/.']) ;
eval(cmd);

cmd =(['! ln -s  '  BiasDir '/Gains.nii.gz ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

%% IV. Organize the T1 weighted images in a directory named 'T1w'

T1wDir=fullfile(mrQ.OutPutNiiDir,'T1w');

% If we redo it and maps directory already exists,
%   we will save the old one before we make a new one.
if (exist(T1wDir,'dir'))
    ex=0; num=0;
    while ex==0
        num=num+1;
        T1wDirOld=fullfile(mrQ.OutPutNiiDir,['T1wDirOld_' num2str(num)] );
        if (~exist(T1wDirOld,'dir'))

            eval(['! mv ' T1wDir  ' ' BiasDirOld]);
            ex=1;
        end
    end
end

mkdir(T1wDir);

%
% cmd =(['! mv ' mrQ.spgr_initDir '/T1wfs_4.nii.gz ' T1wDir '/.']) ;
% eval(cmd);
% 
% cmd =(['! ln -s  '  T1wDir '/T1wfs_4.nii.gz ' mrQ.spgr_initDir '/.']) ;
% eval(cmd);
% 
% cmd =(['! mv ' mrQ.spgr_initDir '/T1wfs_2.nii.gz ' T1wDir '/.']) ;
% eval(cmd);
% 
% cmd =(['! ln -s  '  T1wDir '/T1wfs_2.nii.gz ' mrQ.spgr_initDir '/.']) ;
% eval(cmd);

cmd =(['! mv ' mrQ.T1w_file ' ' T1wDir '/.']) ;
eval(cmd);

cmd =(['! ln -s  '  T1wDir '/T1w.nii.gz ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

cmd =(['! mv ' mrQ.T1w_file1 ' ' T1wDir  '/.']) ;
eval(cmd);

cmd =(['! ln -s  '  T1wDir '/T1w1.nii.gz ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

mrQ.T1w_files=T1wDir;

%% V. Save
 save(mrQ.name,'mrQ');