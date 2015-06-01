function mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir,B1file)
% mrQ_run_Ver2(dir,outDir,useSUNGRID,refFile,inputData_spgr,inputData_seir)
%mrQ_runNIMS(dir,Callproclus,refFile,outDir)
%dir - where the nifti from NIMS are.
% % % % %Callproclus use 1 when using proclus (stanfrod computing cluster)
%refFile a path to a reference image (nii.gz)
%outDir  where the output will be saved (defult is: pwd/mrQ)

%


%%
% Create the initial structure
if notDefined('outDir')
    outDir = fullfile(dir,'mrQ');
end %creates the name of the output directory

if ~exist(outDir,'dir'); mkdir(outDir); end

mrQ = mrQ_Create(dir,[],outDir); %creates the mrQ structure

% Set other parameters
%            mrQ = mrQ_Set(mrQ,'sub',num2str(ii));

if notDefined('useSUNGRID')
    mrQ = mrQ_Set(mrQ,'sungrid',false);
else
    mrQ = mrQ_Set(mrQ,'sungrid',useSUNGRID);
end

%             mrQ = mrQ_Set(mrQ,'sungrid',1);
mrQ = mrQ_Set(mrQ,'fieldstrength',3);

if ~notDefined('refFile')
    mrQ = mrQ_Set(mrQ,'ref',refFile);
    
else
    % New input to automatically acpc align
    mrQ = mrQ_Set(mrQ,'autoacpc',1);
end

%% arrange data
% Specific arrange function for nimsfs nifti or using input for user
if ~isfield(mrQ,'Arange_Date');
    
    if (~notDefined('inputData_spgr') &&  ~notDefined('inputData_seir'))
        mrQ = mrQ_arrangeData_nimsfs(mrQ,inputData_spgr,inputData_seir);
    else
        mrQ = mrQ_arrangeData_nimsfs(mrQ);
        
    end
else
    fprintf('data was allready arrange at %s \n',mrQ.Arange_Date)
end
%%
if notDefined('B1file') 
    % check if B1 was defined by the user. if not we will use the SEIR data
    % to map it.
    
    %% fit SEIR    
    
    if isfield(mrQ,'SEIR_done');
    else
        mrQ.SEIR_done=0;
    end
    
    if (mrQ.SEIR_done==0);
        
        %keep track of the variable we use  for detail see inside the function
        [~, ~, ~, mrQ.SEIRsaveData]=mrQ_initSEIR_ver2(mrQ,mrQ.SEIRepiDir,mrQ.alignFlag);
        
        [mrQ]=mrQ_fitSEIR_T1(mrQ.SEIRepiDir,[],0,mrQ);
        mrQ.SEIR_done=1;
        save(mrQ.name,'mrQ');
        fprintf('fit SEIR  - done!');
        
    else
        fprintf('\n  load fit SEIR data ! \n');
        
    end
    
end
    
    %%
    
