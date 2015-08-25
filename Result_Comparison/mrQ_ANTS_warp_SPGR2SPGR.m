function WarpFiles=mrQ_ANTS_warp_SPGR2SPGR(T1wfile1_LR,T1wfile2_HR,outPutDir,morefiles2_HR)
%
 %mrQ_ANTS_warp_SPGR2SPGR(T1wfile1_LR,T1wfile2_HR,outPutDir,morefiles2_HR)
%inputs:
% T1wfile1_LR          The first T1 weighted image (should be the lower resolution of the two)
 %T1wfile2_HR         The second T1 weighted image  (should be the higher resolution of the two)
% outPutDir             Where to save the registration files and warp parameters
%morefiles2_HR      Other files to be warped 
%% we first need tp make matlab work with the bachrc we have this is
% unfortunate complexity

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

 %% run ANTS non linear parameter to register T1w1 to T1w2 

%name the saved parameter file
 disp('ATNS  linear registration and Warp may take about 20 min ...')
 disp('...')
 
 warpFile= fullfile(outPutDir,'T1w2_to_T1w1');

    % ANTS: make the Ants comand
    cmANTS=['xterm -e ANTS 3 -m CC[' T1wfile1_LR ',' T1wfile2_HR ',1,2] -o ' warpFile '.nii.gz --rigid-affine true'];

     % Run the command in unix and get back status and results: 
    [status, result] = system(cmANTS);
    if status ~= 0 ; 
        cmANTS=['ANTS 3 -m CC[' T1wfile1_LR ',' T1wfile2_HR ',1,2] -o ' warpFile '.nii.gz --rigid-affine true'];
        [status, result] = system(cmANTS,'-echo')
    end


       
%%    name and saved sprg T1w2 map in T1w1 space

% the T1 file name
savefie=fullfile(outPutDir,'Warp_T1w2_to_T1w1.nii.gz');

    % ANTS: make the Ants comand
    cmWarp=['xterm -e WarpImageMultiTransform  3 ' T1wfile2_HR ' ' savefie ' -R '  T1wfile1_LR  ' ' warpFile 'Warp.nii.gz ' warpFile 'Affine.txt'];

    % Run the command in unix and get back status and results:
    [status, result] = system(cmWarp);
    if status ~=0
        cmWarp=['WarpImageMultiTransform  3 ' T1wfile2_HR ' ' savefie ' -R '  T1wfile1_LR  ' ' warpFile 'Warp.nii.gz ' warpFile 'Affine.txt']
        [status, result] = system(cmWarp,'-echo')
    end
        WarpFiles{1}=savefie;
           %%    name and saved raw spgr flipAngles images in epi space  
           
           for d=1:length(morefiles2_HR) %loof over  flip Angles raw images
             file=dir(morefiles2_HR{d});
savefileN=fullfile(outPutDir,['Warp' file.name]);
           
               %ANTS: make the Ants comand
               cmWarp=['xterm -e WarpImageMultiTransform  3 ' morefiles2_HR{d}  ' ' savefileN ' -R '  T1wfile1_LR ' ' warpFile 'Warp.nii.gz ' warpFile 'Affine.txt'];
               % Run the command in unix and get back status and results:
               [status, result] = system(cmWarp);
               if status ~= 0 ;
                   cmWarp=['WarpImageMultiTransform  3 ' morefiles2_HR{d}  ' ' savefileN ' -R '  T1wfile1_LR ' ' warpFile 'Warp.nii.gz ' warpFile 'Affine.txt'];
                   % Run the command in unix and get back status and results:
                   [status, result] = system(cmWarp,'-echo')
               end
               
           WarpFiles{1+d}=savefileN;
            
           end
        
  