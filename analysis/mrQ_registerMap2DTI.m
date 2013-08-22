function mrQ_registerMap2DTI(B0file,T1file,otherMaps,outDir)


orig_path = getenv('LD_LIBRARY_PATH');

colon_idx = strfind(orig_path, ':');

matlab_path = [];
other_path = [];

for i=1:(length(colon_idx)-1)
    
    this_path = orig_path(colon_idx(i):colon_idx(i+1));
    % This is part of the matlab path:
    if strfind(this_path, 'matlab')
        matlab_path = [matlab_path, ':', this_path];
        
        % Should go before the matlab bit of the path (all the rest):
    else
        other_path = [other_path, ':', this_path];
        
    end
end

new_path = [other_path, ':', matlab_path];

setenv('LD_LIBRARY_PATH', new_path)

out=fullfile(outDir,'map2B0');


%make a trasformation from T1file (moving image) to Bofile (fix image)

cmANTS=['xterm -e ANTS 3 -m CC[' B0file ',' T1file ',1,2] -o ' out '.nii.gz --rigid-affine true'];

% Run the command in unix and get back status and results:
[status, result] = system(cmANTS);

if status ~= 0

    cmANTS=['ANTS 3 -m CC[' B0file ',' T1file ',1,2] -o ' out '.nii.gz --rigid-affine true'];
    
    % Run the command in unix and get back status and results:
    [status, result] = system(cmANTS,'-echo');
end



%applay the trasformation to the images

[~, name]=fileparts(T1file);
[~, name]=fileparts(name);

out1=fullfile(outDir,[name '_2DTI.nii.gz']);

cmWarp=['xterm -e WarpImageMultiTransform  3 ' T1file  ' ' out1 ' -R '  T1file ' ' out 'Warp.nii.gz ' out 'Affine.txt' ' --use-BSpline'];

% Run the command in unix and get back status and results:
[status, result] = system(cmWarp);
% register the other maps
if ~notDefined('otherMaps')
    for i=1:length(otherMaps)
        [~, name]=fileparts(otherMaps{i});
        [~, name]=fileparts(name);
        out1=fullfile(outDir,[name '_2DTI.nii.gz']);
        
        cmWarp=['xterm -e WarpImageMultiTransform  3 ' otherMaps{i}  ' ' out1 ' -R '  T1file ' ' out 'Warp.nii.gz ' out 'Affine.txt' ' --use-BSpline'];
        % Run the command in unix and get back status and results:
        
        [status, result] = system(cmWarp);
        
        if status ~= 0
            cmWarp=['WarpImageMultiTransform  3 ' otherMaps{i}  ' ' out1 ' -R '  T1file ' ' out 'Warp.nii.gz ' out 'Affine.txt' ' --use-BSpline'];
            % Run the command in unix and get back status and results:
            
            [status, result] = system(cmWarp,'-echo');
        end
    end
end



%%
% get the path back

setenv('LD_LIBRARY_PATH', orig_path)
