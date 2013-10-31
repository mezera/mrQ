function SIEMENS_dcm2nii(dicomDir, combineddicomDir, savefilename)

% Preliminary script to sort T1 flash dicom file in SIEMENS scanner in
% order to create a nifti file used for mrQ analyses.
%
% INPUT:
% dicomDir: a directory including DICOM for t1 flash scan in each channel
% combineddicomDir: a directory including DICOM combining t1 flash scan in all channels
% savefilename: output nifti file name

%
% Hiromasa Takemura 2013 (C) Stanford VISTA team

%% Read dicom file names
d = genpath(dicomDir);

% Argument checking about dicom direcctory
if(isempty(d))
    error(['Dicom dir "' d '" not found or empty.']);
end
if(isunix)
    d = explode(':',genpath(dicomDir));
else
    d = explode(';',genpath(dicomDir));
end
d = d(1:end-1);
n = 0;

% Specifying the file name of DICOMs
for(ii=1:length(d))
    d2 = dir(fullfile(d{ii},'*'));
    for(jj=1:length(d2))
        if(~d2(jj).isdir)
            n = n+1;
            dicomFiles{n} = fullfile(d{ii},d2(jj).name);
        end
    end
end

clear d;
numSeries = 0;

%% Read the dicom info
for(imNum=1:length(dicomFiles))
    curFile = dicomFiles{imNum};
    
    % Read Dicom header in each file
    info = dicominfo(curFile);
    
    
    
    
    
    
    
    %%
    % TR
    s.TR = info.RepetitionTime;
    s.seriesDescription = info.SeriesDescription;
    
    % Phase encoding direction
    s.phaseEncodeDir = info.InPlanePhaseEncodingDirection;
    
    % Image orientation
    s.imageOrientation = info.ImageOrientationPatient;
    
    % Spatial resolution of data
    if(isfield(info,'SpacingBetweenSlices'))
        s.mmPerVox = [info.PixelSpacing(:)' info.SpacingBetweenSlices];
    else
        s.mmPerVox = [info.PixelSpacing(:)' info.SliceThickness];
    end
    
    
    % Number of slices
    sliceNum = info.InstanceNumber;
    
    % Image Position
    s.imagePosition(:,sliceNum) = info.ImagePositionPatient;
    
    % Slice Location
    s.sliceLoc(sliceNum) = info.SliceLocation;
    %flip angle
    s.flipAngle = info.FlipAngle;
    % TE
    s.TE=info.EchoTime;


    if (isfield(info,'InversionTime'))
        s.inversionTime=info.InversionTime;
    else
        s.inversionTime=[];
    end
        


    % Read the header which specifies which dicom is from which channel
    channeltag = info.Private_0051_100f;
    [token,remain] = strtok(channeltag,'H');
    channelcode = str2num(token);
    s.imData(:,:,sliceNum, channelcode) = dicomread(info)';
end

%% Read the dicom files combining all slices

% Argument chekinc of dicom folder
d = genpath(combineddicomDir);
if(isempty(d))
    error(['Dicom dir "' d '" not found or empty.']);
end
if(isunix)
    d = explode(':',genpath(combineddicomDir));
else
    d = explode(';',genpath(combineddicomDir));
end
d = d(1:end-1);
n = 0;

% Specifying the dicom file name
for(ii=1:length(d))
    d2 = dir(fullfile(d{ii},'*'));
    for(jj=1:length(d2))
        if(~d2(jj).isdir)
            n = n+1;
            dicomFiles2{n} = fullfile(d{ii},d2(jj).name);
        end
    end
end

clear d;
numSeries = 0;

% Attaching combined dicom info to nifti file as a 33th channel
for(imNum=1:length(dicomFiles2))
    curFile = dicomFiles2{imNum};
    info = dicominfo(curFile);
    s.TR = info.RepetitionTime;
    s.seriesDescription = info.SeriesDescription;
    s.phaseEncodeDir = info.InPlanePhaseEncodingDirection;
    sliceNum = info.InstanceNumber;
    s.imData(:,:,sliceNum, 33) = dicomread(info)';
end

%% Compute scanner-to-image xform

dim = size(s.imData);

% Orientation information
%-------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system:
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% T&T co-ordinate system:
% x increases      left to right
% y increases posterior to anterior
% z increases  inferior to superior

% Code lifted from dinifti source:
% patient orientation cosines, in [Sag, Cor, Ax] order
rowCosines = s.imageOrientation(1:3);
colCosines = s.imageOrientation(4:6);
sliceNorm = cross(rowCosines, colCosines);

qto_xyz = zeros(4,4);
qto_xyz(1,1) = -rowCosines(1); % -image->RowSagCos();
qto_xyz(1,2) = -colCosines(1); % -image->ColSagCos();
qto_xyz(1,3) = -sliceNorm(1); % -(rowCosines(2) * colCosines(3) - rowCosines(3) * colCosines(2)); % -image->SagNorm();

qto_xyz(2,1) = -rowCosines(2); % -image->RowCorCos();
qto_xyz(2,2) = -colCosines(2); % -image->ColCorCos();
qto_xyz(2,3) = -sliceNorm(2); % -(rowCosines(3) * colCosines(1) - rowCosines(1) * colCosines(3)); % -image->CorNorm();

qto_xyz(3,1) = rowCosines(3); % image->RowTraCos();
qto_xyz(3,2) = colCosines(3); % image->ColTraCos();
qto_xyz(3,3) = sliceNorm(3); % rowCosines(1) * colCosines(2) - rowCosines(2) * colCosines(1); % image->TraNorm();


% Simple hack to decide if we need to flip the slice order:
if(dot(sliceNorm,s.imagePosition(:,1)) > dot(sliceNorm,s.imagePosition(:,end)))
    %         s.sliceNum = flipdim(s.sliceNum,2);
    s.sliceLoc = flipdim(s.sliceLoc,2);
    s.imagePosition = flipdim(s.imagePosition,2);
    s.imData = flipdim(s.imData,3);
end

% offset = slicePosition(center) - M*(i,j,k of center) since we are working with the first slice k = 0.

% For matlab's 1-indexing, we have to subtract mmPerVox from the origin
% (imagePosition) to reflect the fact that the first voxel is [1,1,1]
% rather than [0,0,0].
pos = s.imagePosition(:,1);
qto_xyz(:,4) = [-pos(1) -pos(2) pos(3) 1]';
qto_xyz(1:3,1:3) = qto_xyz(1:3,1:3)*diag(s.mmPerVox);

% Apply a 1-voxel offset to the origin to account for Matlab's
% 1-indexing.
qto_xyz = inv(inv(qto_xyz)+[0 0 0 1;0 0 0 1;0 0 0 1; 0 0 0 0]);


s.imToScanXform  = qto_xyz;


%% Write nifti
if(strcmpi(s.phaseEncodeDir,'ROW')) fpsDim = [2 1 3];
else fpsDim = [1 2 3]; end

% Initilize nifti structure
ni = niftiGetStruct(s.imData, s.imToScanXform, 1, s.seriesDescription, [], [], fpsDim);

%for now we are keeping the same trasform to the nifti and dicoms
% Applying xform to nifti  
ni = niftiApplyCannonicalXform(ni);

% Setting nifti file name
ni.fname = savefilename;

% Save file
niiName=[savefilename '.nii.gz'];
niftiWriteMatlab(ni,savefilename);
strName=[savefilename '.mat'];

save(strName,'s')

