function mrQ=mrQ_arangeOutPutDir(mrQ)
fprintf('\n orgenized the output directory                \n');

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


%%  orgnized the brain maps in a directory named maps
fprintf('\n orgenized the maps in a maps directory                \n');
mapDir=fullfile(mrQ.OutPutNiiDir,'BrainMaps');

% if we redo it and maps dir is already exsist we will saved
% the old before we make a new one
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
  T1file=mrQ_getT1file(mrQ);
  
  [d dd ddd]=fileparts(T1file);
 cmd =(['! mv ' T1file ' ' mapDir '/.']);
     eval(cmd);
  
    cmd =(['! ln -s  '  mapDir '/' dd ddd ' ' mrQ.spgr_initDir '/.']) ;
eval(cmd);
% % If T1_map_lsq.nii.gz exists in the spgr_intDir then move it to the map
% % dir. Otherwise assume that it ended up in mapDirOld and retain this t1
% % map. This would occur if the pd fits were rerun but the t1 fits were not
% if exist(fullfile(mrQ.spgr_initDir, 'T1_map_lsq.nii.gz'),'file')
%     cmd =(['! mv ' mrQ.spgr_initDir '/T1_map_lsq.nii.gz ' mapDir '/.']);
%     eval(cmd);
%     
%     cmd =(['! ln -s  '  mapDir '/T1_map_lsq.nii.gz ' mrQ.spgr_initDir '/.']) ;
% eval(cmd);
% 
% else
%     % If a new t1 map was not generated than copy it from the old
%     % directory. We copy rather than move to leave the old directory intact
%     cmd =(['! cp ' mapDirOld '/T1_map_lsq.nii.gz ' mapDir '/.']);
%     eval(cmd);
% end


cmd =(['! mv ' mrQ.spgr_initDir '/WF_map.nii.gz ' mapDir '/.']) ;
eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/VIP_map.nii.gz ' mapDir '/.']) ;
eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/TV_map.nii.gz ' mapDir '/.']) ;
eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/SIR_map.nii.gz ' mapDir '/.']) ;
eval(cmd);



mrQ.mapsDir=mapDir;
mrQ.maps.T1path=fullfile(mapDir,'T1_map_lsq.nii.gz');
mrQ.maps.WFpath=fullfile(mapDir,'WF_map.nii.gz');
mrQ.maps.TVpath=fullfile(mapDir,'TV_map.nii.gz');
mrQ.maps.SIRpath=fullfile(mapDir,'SIR_map.nii.gz');

%% SEIR

SEIRT1Dir=fullfile(mapDir,'SEIR');
mkdir(SEIRT1Dir)
cmd =(['! cp ' mrQ.T1file ' '  SEIRT1Dir '/.']) ;
eval(cmd);



%% BIAS

BiasDir=fullfile(mrQ.OutPutNiiDir,'BiasMap');

% if we redo it and maps dir is already exsist we will saved
% the old before we make a new one
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

%% wighted images

T1wDir=fullfile(mrQ.OutPutNiiDir,'T1w');

% if we redo it and maps dir is already exsist we will saved
% the old before we make a new one
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
cmd =(['! mv ' mrQ.spgr_initDir '/T1wfs_4.nii.gz ' T1wDir '/.']) ;
eval(cmd);

cmd =(['! ln -s  '  T1wDir '/T1wfs_4.nii.gz ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

cmd =(['! mv ' mrQ.spgr_initDir '/T1wfs_2.nii.gz ' T1wDir '/.']) ;
eval(cmd);

cmd =(['! ln -s  '  T1wDir '/T1wfs_2.nii.gz ' mrQ.spgr_initDir '/.']) ;
eval(cmd);

%%
 save(mrQ.name,'mrQ');