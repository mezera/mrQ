function alignedMaps = mrQ_registerMap2DTI(B0file,T1file,otherMaps,outDir, interpMethod, resampleFlag,recalc)
%
%
% alignedMaps = mrQ_registerMap2DTI(B0file,T1file,otherMaps,outDir, interpMethod)
%
%
% resampleFla is a binary flag:
%   - resampleFlag = 0 (default) will result in T1 (and other maps) staying in
%     their native resolution.
%   - resampleFlag = 1 will result in resampling of the T1file (and other maps)
%     into the diffusion space resolution.
%
% (C) Stanford University, VISTA LAB - 2014


%% Handle INPUTS
if notDefined('interpMethod')
    interpMethod = '--use-BSpline';
elseif strcmp(interpMethod,'spline')
    interpMethod = '--use-BSpline';
elseif strcmp(interpMethod,'nn') || strcmp(interpMethod,'nearestneighbor')
    interpMethod =  '--use-NN';
end

if notDefined('resampleFlag')
    % by default, no resampling to the reference map is done (T1 will stay
    % in its originla resolution)
    resampleFlag = 0;
end
if notDefined('recalc')
    recalc=1;
end

% In case of no resampling, the transformed maps will have the same
% resolution as the T1File. Otherwise, they will have same resolution as
% the B0File.
if resampleFlag == 0
    refFile = T1file;
else
    refFile = B0file;
end


%% Configure ENV paths
% %
% for i=1:(length(colon_idx)-1)
%
%     this_path = orig_path(colon_idx(i):colon_idx(i+1));
%     % This is part of the matlab path:
%     if strfind(this_path, 'matlab')
%         matlab_path = [matlab_path, ':', this_path];
%
%         % Should go before the matlab bit of the path (all the rest):
%     else
%         other_path = [other_path, ':', this_path];
%
%     end
% end
%
% new_path = [other_path, ':', matlab_path];
%
% setenv('LD_LIBRARY_PATH', new_path)
%
%

%% Calculate the transform

% Set the name for the output of the function
out=fullfile(outDir,'map2B0');

% make a trasformation from T1file (moving image) to Bofile (fixed image)
cmANTS=['xterm -e ANTS 3 -m CC[' B0file ',' T1file ',1,2] -o ' out '.nii.gz --rigid-affine true'];

if recalc
% Run the command in unix and get back status and results:
[status, ~] = system(cmANTS);

% If the above command fails...
if status ~= 0
    cmANTS=['ANTS 3 -m CC[' B0file ',' T1file ',1,2] -o ' out '.nii.gz --rigid-affine true'];
    % Run the command in unix and get back status and results:
    [~, ~] = system(cmANTS,'-echo');
end
end



%% Apply the trasformation to the images
[~, name] = fileparts(T1file);
[~, name] = fileparts(name);

if resampleFlag == 0
    out1 = fullfile(outDir,[name '_2DTI.nii.gz']);
else
    out1 = fullfile(outDir,[name '_2DTI_resamp.nii.gz']);
end


 cmWarp=['xterm -e WarpImageMultiTransform  3 ' T1file  ' ' out1 ' -R ' refFile ' ' out 'Warp.nii.gz ' out 'Affine.txt ' interpMethod];
% Run the command in unix and get back status and results:
 [status, ~] = system(cmWarp);

 if status ~= 0 
      cmWarp=['WarpImageMultiTransform  3 ' T1file  ' ' out1 ' -R ' refFile ' ' out 'Warp.nii.gz ' out 'Affine.txt ' interpMethod]; 
      [~, ~] = system(cmWarp,'-echo');
  end
% Start returning the aligned files
alignedMaps{1} = out1;


%% Register the other maps

if ~notDefined('otherMaps')
    for i=1:length(otherMaps)
        [~, name]=fileparts(otherMaps{i});
        [~, name]=fileparts(name);
        
        if resampleFlag == 0
            out1=fullfile(outDir,[name '_2DTI.nii.gz']);
        else
            out1=fullfile(outDir,[name '_2DTI_resamp.nii.gz']);
        end
        
        cmWarp=['xterm -e WarpImageMultiTransform  3 ' otherMaps{i}  ' ' out1 ' -R ' refFile ' ' out 'Warp.nii.gz ' out 'Affine.txt ' interpMethod];
        
        % Run the command in unix and get back status and results:
        [status, ~] = system(cmWarp);
        
        if status ~= 0
            cmWarp=['WarpImageMultiTransform  3 ' otherMaps{i}  ' ' out1 ' -R ' refFile ' ' out 'Warp.nii.gz ' out 'Affine.txt ' interpMethod];
            % Run the command in unix and get back status and results:
            [~, ~] = system(cmWarp,'-echo');
        end
        
        % Return the aligned maps
        alignedMaps{end+1} = out1;
    end
end



%%

% get the path back
% setenv('LD_LIBRARY_PATH', orig_path)

return


