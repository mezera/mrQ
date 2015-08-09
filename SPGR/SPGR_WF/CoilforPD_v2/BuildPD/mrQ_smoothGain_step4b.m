function [opt, G, PD]=mrQ_smoothGain_step4b(opt,PD,PD_filename,G_filename)
% [opt, G, PD]=mrQ_smoothGain_step4b(opt,PD,PD_filename,G_filename)
%
% This is Step 5 of 6 (including Step 0) in the pipeline to build the WF
% (water fraction) map. In this step, the function gets a smooth coil
% sensitivity in all locations, brings it back to the original image space,
% and calculates the PD.
%
%
% ~INPUTS~
%             opt:   The "opt" structure of optimized parameters
%              PD:
%     PD_filename:   The PD NIfTI file as stored in the opt structure. If 
%                          no PD_filename is selected, the default is to 
%                          take from opt's outDir, PD.nii.gz
%      G_filename:   The Gains NIfTI file as stored in the opt structure. 
%                          If no G_filename is selected, the default is to 
%                          take from opt's outDir, Gains.nii.gz
%
% ~OUTPUTS~
%             opt:   The updated opt structure of optimized parameters
%               G:   Location of the newly constructed Gain image
%              PD:   Location of the newly constructed PD image
%
%
% See also: mrQ_buildPD_ver2
%           Step_0: none
%           Step_1: mrQ_CalBoxPD_step1a
%           Step_2: mrQ_ScaleBoxes_step2
%           Step_3: mrQ_BoxJoinBox
%           Step_5: mrQ_PD2WF_step5
%
%
% (C) Mezer lab, the Hebrew University of Jerusalem, Israel
%   2015
%
%


%% I. Load NIfTI files and filenames

%multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;

%M0=sqrt(sum(M0.^2,4));
M0=(sum(M0,4));

%Brain mask
BM=readFileNifti(opt.BMfile);
xform=BM.qto_xyz;
BM=BM.data;

if notDefined('PD_filename');PD_filename=fullfile(opt.outDir, 'PD.nii.gz');end
if notDefined('G_filename');G_filename=fullfile(opt.outDir, 'Gains.nii.gz');end

%% II. Calculate coil gains from the PD_fit and M0
% M0 is a 4D image, where the 4th dimension is the image of each coil:
%        M0= Gain * PD  --> Gain= M0 / PD;
%
% PD fit is not a full image: We have values only where the box fits were
% done and a "good solution" was obtained. (Hopefully we have most of the
% image.) 
% It is given that the Gain function varies smoothly in space. Therefore,
% to fill the areas where PD data is missing, we will model the Gain as a
% smooth function in space.
%
% This model will also smoothe areas where the Gain estimation is noisy.
mask=BM & PD>0 ;
PD=PD./median(PD(mask));
Gi=zeros(size(PD));
Gi=M0./PD;
Gi(~mask)=nan;
Gi(isinf(Gi))=nan;

%% III-a. Smoothe by looping over the Z slices

% We will estimate a smooth function in 2D space (many Z-slices)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,2));

Gz=zeros(size(M0));

%run over coils and fit the Gain for each

%loop over slices
for  jj=1:size(PD,3)
    
    tmp=Gi(:,:,jj);
    tmpBM=BM(:,:,jj);
    
    %Check that there is data in the slice
    wh=find(tmp>0);
    if    (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>800)
        
        %Find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
        ZI = griddata(xg,yg,zg,XI,YI);
        
        % Put the resulting Gain in the 3D Gain image and fix the orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        Gz(:,:,jj)=ZI;
        
        clear ZI
    end;  
    
end

PDz=M0./Gz;
PDz(PDz>2)=0; PDz(PDz<0)=0;

%% III-b. Smoothe by looping over the Y slices

% We will estimate a smooth function in 2D space (many Y-slices)

[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,3));

Gy=zeros(size(M0));

%run over coils and fit the Gain for each

%loop over slices
for  jj=1:size(PD,2)
    
    tmp=squeeze(Gi(:,jj,:));
    tmpBM=squeeze(BM(:,jj,:));
    
    %Check that there is data in the slice
    wh=find(tmp>0);
    if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>800);
        
        %Find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
        ZI = griddata(xg,yg,zg,XI,YI);
        
        % Put the resulting Gain in the 3D Gain image and fix the orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        Gy(:,jj,:)=ZI;
        
        clear ZI
    end;
    
end;

PDy=M0./Gy;
PDy(PDy>2)=0; PDy(PDy<0)=0;

%% III-c. Smoothe by looping over the X slices

% We will estimate a smooth function in 2D space (many X-slices)
[XI, YI]=meshgrid(1:size(PD,2),1:size(PD,3));

Gx=zeros(size(M0));

%run over coils and fit the Gain for each

%loop over slices
for  jj=1:size(PD,1)
    
    tmp=squeeze(Gi(jj,:,:));
    tmpBM=squeeze(BM(jj,:,:));
    
    %Check that there is data in the slice
    wh=find(tmp>0);
    if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>800);
        
        %Find location of data
        [x,y] = ind2sub(size(tmp),wh);
        z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
        
        [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
        ZI = griddata(xg,yg,zg,XI,YI);
        
        % Put the resulting Gain in the 3D Gain image and fix the orientation
        ZI=rot90(ZI);
        ZI = flipdim(ZI,1);
        Gx(jj,:,:)=ZI;
        
        clear ZI
    end;
    
end;

PDx=M0./Gx;
PDx(PDx>2)=0; PDx(PDx<0)=0;

%% IV. Join the three dimensions' PD, and calculate the Gain

PDfit=zeros(size(PD));
PDvox=zeros(size(PD));

whx=find(PDx>0); %x
PDfit(whx)=PDx(whx);
PDvox(whx)=PDvox(whx)+1;


why=find(PDy>0); %y
PDfit(why)=PDfit(why)+PDy(why);
PDvox(why)=PDvox(why)+1;



whz=find(PDz>0); %z
PDfit(whz)=PDfit(whz)+PDz(whz);
PDvox(whz)=PDvox(whz)+1;

PDfit=PDfit./PDvox;

G=M0./PDfit;

G(isinf(G))=nan;
G(G<0)=nan;

Gfit=zeros(size(M0));
mask=PDfit>0 & BM;

%% V-a. Estimate Gz
% We will estimate a smooth function in 2D space (many Z-slices)
[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,2));

%loop over slices
    for  jj=1:size(PD,3)
        
        tmp=G(:,:,jj);
        tmpBM=mask(:,:,jj);
        
        %Check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(find(tmpBM>0))>0.3  && length(wh)>1000);
            
            %Find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            
        % Put the resulting Gain in the 3D Gain image and fix the orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gfit(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
%% V-b. Estimate Gx
% We will estimate a smooth function in 2D space (many X-slices)

[XI, YI]=meshgrid(1:size(PD,2),1:size(PD,3));

%loop over slices
    for  jj=1:size(PD,1)
        
        tmp=squeeze(Gfit(jj,:,:));
        
        %Check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmpBM(:))<1  && length(wh)>800);
            %Find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            
        % Put the resulting Gain in the 3D Gain image and fix the orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gfit(jj,:,:)=ZI;
            
            clear ZI
        end;
        
    end;

%% V-c. Estimate Gy
% We will estimate a smooth function in 2D space (many Y-slices)

[XI, YI]=meshgrid(1:size(PD,1),1:size(PD,3));

%loop over slices
    for  jj=1:size(PD,2)
        
        tmp=squeeze(Gfit(:,jj,:));
        
        %Check that there is data in the slice
        wh=find(tmp>0);
        if   (length(find(tmp>0))/length(tmp(:))>0.3 && length(find(tmp>0))/length(tmpBM(:))<1  && length(wh)>800);
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
        % Estimate a smooth version of the data in the slice. 
        % For original code see: 
        %           Noterdaeme et al. "Intensity correction with a pair of 
        %               spoiled gradient recalled echo images". Phys. Med. 
        %               Biol. 54 3473-3489 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            
        % Put the resulting Gain in the 3D Gain image and fix the orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gfit(:,jj,:)=ZI;
            
            clear ZI
        end;
        
    end;
   
Gfit(isinf(Gfit))=0;
Gfit(Gfit<0)=0;
Gfit(isnan(Gfit))=0;

%% VI. Moving from the Gain back to the PD
PDfit1=M0./Gfit;

PDfit1(isinf(PDfit1))=0;
PDfit1(PDfit1<0)=0;
PDfit1(isnan(PDfit1))=0;
PDfit1(PDfit1>4)=4;

%% VII. Save the PD and Coil gain
% Upsample if needed to the original resolution

% Save the Gain in the resolution for which we calculated it
dtiWriteNiftiWrapper(single(Gfit),xform,G_filename);

if ~isfield(opt,'Resamp')
    opt.Resamp=0;
end
if ~opt.Resamp==1; %no resample; we will save PD and coil Gain
    dtiWriteNiftiWrapper(single(PDfit1),xform,PD_filename);
    
else
    % We will upsample the Gain to the original M0 size and then divide
    %   (e.g., PD=M0/G)
    
    sz=size(Gfit);
    bb = mrAnatXformCoords(xform, [1 1 1; sz(1:3)]);
    
    % get the original M0
    clear M0
    M0=readFileNifti(opt.M0file_Org);outMm=M0.pixdim(1:3);     M0=M0.data; MSZ=size(M0);
    M0=(sum(M0,4));
    
    PD=zeros(size(M0));
    [Gi, M0UnderSamp_Xform] = mrAnatResliceSpm(double(Gfit(:,:,:)), inv(xform), bb, outMm, 1, 0);
    GSZ=size(Gi);
    % Get the size difference, if it exists
    Sizes=[min(GSZ(1),MSZ(1)) min(GSZ(2),MSZ(2)) min(GSZ(3),MSZ(3)) ]; % we calculate the sizes of the two image M0 and G becouse resample may end up with adding or losing a line of zeros.
    
    % Calculate the coil PD
    PD(1:Sizes(1),1:Sizes(2),1:Sizes(3))=M0(1:Sizes(1),1:Sizes(2),1:Sizes(3))./Gi(1:Sizes(1),1:Sizes(2),1:Sizes(3));
    PD(isinf(PD))=0;
    PD(PD<0)=0;
    PD(isnan(PD))=0;
    PD(PD>4)=4;
    
    Gii=zeros(MSZ(1:3));
    Gii(1:Sizes(1),1:Sizes(2),1:Sizes(3))=Gi(1:Sizes(1),1:Sizes(2),1:Sizes(3));
    
    dtiWriteNiftiWrapper(single(Gii),M0UnderSamp_Xform,G_filename);
    
    dtiWriteNiftiWrapper(single(PD),M0UnderSamp_Xform,PD_filename);
    
end

opt.PDfile=PD_filename;
opt.Gainfile=G_filename;

%save(opt.logname,'opt')
