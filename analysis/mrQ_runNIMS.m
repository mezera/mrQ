function mrQ_runNIMS(dir,Callproclus,refFile,outDir)


   % Create the initial structure
   if notDefined('outDir')
            outDir = fullfile(dir,'mrQ');
   end
            if ~exist(outDir,'dir'); mkdir(outDir); end
            mrQ = mrQ_Create(dir,[],outDir);

            % Set other parameters
%            mrQ = mrQ_Set(mrQ,'sub',num2str(ii));
            
            if notDefined('Callproclus')
                mrQ = mrQ_Set(mrQ,'proclus',false);
            else
                mrQ = mrQ_Set(mrQ,'proclus',Callproclus);
            end
            
            mrQ = mrQ_Set(mrQ,'sungrid',1);
            mrQ = mrQ_Set(mrQ,'fieldstrength',3);

      

            % Specific arrange function for nimsfs
            mrQ = mrQ_arrangeData_nimsfs(mrQ);
            
            if ~notDefined('refFile')
                mrQ = mrQ_Set(mrQ,'ref',refFile);
                
            else
                % New input to automatically acpc align
                mrQ = mrQ_Set(mrQ,'autoacpc',1);
            end
            
            % RUN IT
            mrQ_run(mrQ.name);