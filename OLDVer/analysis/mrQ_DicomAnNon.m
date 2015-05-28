function mrQ_DicomAnNon(dicomsDir,delOriginal)
% this function will annonamizde all dicoms in the  dicomsDir
% it will remove the original dicom unless spesicied delOriginal=0


if notDefined('delOriginal')
    delOriginal=true ;
end

%find the dir in th path
dcms = dir(dicomsDir);

% number of directories
NDir=size(dcms,1);
%
for ii=3:NDir
    

        dicomfile=     (fullfile(dicomsDir, dcms(ii).name));
        dicomfileAN= (fullfile(dicomsDir, ['AN_' dcms(ii).name]));
        dicomanon(dicomfile,dicomfileAN)
        if delOriginal
            eval( ['! rm ' dicomfile])
        end
 
end