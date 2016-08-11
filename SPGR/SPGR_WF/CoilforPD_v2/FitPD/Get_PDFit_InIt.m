function[ PDinit,g0] =Get_PDFit_InIt(Init,M0,pBasis,R1,seg)

% Init           approch to init the search  1 defult (SOS) , 2 T1 PD, theortical relations, 3 segmentation

%[ PDinit,g0] =Get_PDFit_INIt(Ref,M0,pBasis,R1,seg)
% estimate a PD from the data

% thuis function define the way we guess the PD (PDinit) and polyinmyal coefisent (g0) to intiate the search
% Init define the method we will use
%  Init=1 or using sum of sqaure (SOS) (defult)
%  Init=2 using the R1 PD theortical relationship
% Init=3 useing a segmentation image taking assuming that a single tissue type (the most abanded) in the segmentation image have a constant value.

% AM  & BW VISTASOFT Team, 2013

if notDefined('Init')
    Init=1;
end


%%

if  Init==1 % sum of squere
    PDinit = sqrt(sum(M0.^2,2));
end


if  Init==2 %R1 PD linear relations
    
    PDinit=1./(R1*0.42+0.95); %this is the teortical T1 PD relationship see reference at Mezer et. al 2013
end

if Init==3 % use segmentation
    TissueType=unique(seg);
    for ii=1:max(TissueType);
        Nvox(ii)=length(find(seg==ii));
    end
    [~, ind]=sort(Nvox);
    PDinit=zeros(size(seg));
    PDinit(seg==ind(end))=1;
    PDinit=PDinit(:);
end



%guess the g0

if ~notDefined('pBasis')
    
    nVoxels=size(M0,1);
    Ncoils=size(M0,2);
    nPolyCoef=size(pBasis,2);
    
    G  = zeros(nVoxels,Ncoils);
    g0 = zeros(nPolyCoef,Ncoils);
    mask=find(PDinit~=0);
    for ii=1:Ncoils
        G(mask,ii)  = M0(mask,ii) ./ PDinit(mask);         % Raw estimate
        g0(:,ii) = pBasis(mask,:) \ G(mask,ii);  % Polynomial approximation
    end
else
    g0=[];
end