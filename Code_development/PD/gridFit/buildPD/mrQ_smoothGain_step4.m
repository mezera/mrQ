function [opt, G, PD]=mrQ_smoothGain_step4(opt,PD)

%multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;
%T1

%Brain mask
BM=readFileNifti(opt.BMfile);
xform=BM.qto_xyz;
BM=BM.data;

%%
%calculate coil gains from the PD_fit and M0. M0 is 4D image the 4th
%dimention is the image of each coil
% M0=Gain X PD  --> Gain= M0 / PD;
% PD fit is not full image. we have values only where the boxes fit where done
% and got a "good sulotion". hopfully we have most of the image.
% It is given that the gain function vart smoothly in space.
% Therefore to fill the area where PD data is missing we will model the
% Gain as a smooth function in space.


% This model will also smooth area where the gain estimation is noise.
mask=BM & PD>0 ;

% we will estimate  a smooth function in 2D space (Z -slice)
[XI YI]=meshgrid(1:size(PD,1),1:size(PD,2));

G=zeros(size(M0));

% we will estimate any other spcae by global smooth poly in the end
Imsz=size(mask);degree=5;
[Poly1,str] = constructpolynomialmatrix3d(Imsz,find(ones(Imsz)),degree);


%l run over coils and fit the Gain for each
for ii=1:size(M0,4)
    
    
    %the estimate gain for the coil
    Gi=zeros(size(PD));
    Gi=M0(:,:,:,ii)./PD;
    Gi(~mask)=nan;
    Gi(isinf(Gi))=nan;
    
    
    %loop over slices
    for  jj=1:size(PD,3)
        
        tmp=Gi(:,:,jj);
        
        %check that there is data in the slice
        wh=find(tmp>0);
        if  length(find(tmp>0))>100;
            %find location of data
            [x,y] = ind2sub(size(tmp),wh);
            z=double(tmp(wh));
            % estimate a smooth vertion of the data in the slice for original code see:
            % Moterdaeme et.al. Phys. Med. Biol. 54 3474-89 (2009)
            
            [zg,xg,yg]= gridfit(x,y,z,1:2:size(tmp,1),1:2:size(tmp,2),'smoothness',5);
            ZI = griddata(xg,yg,zg,XI,YI);
            % put the result gain in the 3D gain image and fix orientation
            ZI=rot90(ZI);
            ZI = flipdim(ZI,1);
            Gi(:,:,jj)=ZI;
            
            clear ZI
        end;
        
    end;
    
    clear param Gitmp maskG
    %fit the polynials coefitents to the smooth Gi map
    [params] = fit3dpolynomialmodel(Gi,(Gi>0),degree);
    
    %reshape from a vector to 3D map
    Gitmp= reshape(Poly1*params(:),Imsz);
    
    %find where the Gi value are missing
    maskG=logical(Gi>0);
    
    %fill the holls
    Gi(~maskG)=Gitmp(~maskG);
    
    
    %save the Gain result for this coil
    G(:,:,:,ii)=Gi;
    
end


PD=M0./G;
PD=median(PD,4);

%% save PD and Coil gain & Upsample if needed to the original resulotion

% save the gain in the resulotion we calculate it
G_filename=fullfile(opt.outDir, 'Gains.nii.gz');
dtiWriteNiftiWrapper(single(G),xform,G_filename);

PD_filename=fullfile(opt.outDir, 'PD.nii.gz');


if ~isfield(opt,'Resamp')
    opt.Resamp=0;
end
if ~opt.Resamp==1; %no resample we will save PD and coil gain
    dtiWriteNiftiWrapper(single(PD),xform,PD_filename);
    
else
    %we will up sample the G to the original M0 size and then divied M0/G
    %to get PD
    sz=size(G);
    bb = mrAnatXformCoords(xform, [1 1 1; sz(1:3)]);
    clear M0
    % get the original M0
    M0=readFileNifti(opt.M0file_Org);outMm=M0.pixdim(1:3);     M0=M0.data; MSZ=size(M0);
    % loop over coils and up sample G
    for ii=1:sz(4)
        
        [Gi, M0UnderSamp_Xform] = mrAnatResliceSpm(double(G(:,:,:,ii)), inv(xform), bb, outMm, 1, 0);
        GSZ=size(Gi);
        % get the size difference is exsist
        Sizes=[min(GSZ(1),MSZ(1)) min(GSZ(2),MSZ(2)) min(GSZ(3),MSZ(3)) ]; % we calculate the sizes of the two image M0 and G becouse resample may end up with adding or losing a line of zeros.
        
        % calculate the coil PD
        M0(1:Sizes(1),1:Sizes(2),1:Sizes(3),ii)=M0(1:Sizes(1),1:Sizes(2),1:Sizes(3),ii)./Gi(1:Sizes(1),1:Sizes(2),1:Sizes(3));
        
    end
    PD=median(M0,4); %is wight sum by SNR better
    
    dtiWriteNiftiWrapper(single(PD),M0UnderSamp_Xform,PD_filename);
    
end

opt.PDfile=PD_filename;
opt.Gainfile=G_filename;
save(opt.logname,'opt')