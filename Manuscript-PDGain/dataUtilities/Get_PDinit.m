function[ PDinit,g0] =Get_PDinit(T1reg,R1,Guesstype,M0,pBasis)
%[ PDinit,g0] =Get_PDinit(T1reg,R1,Guesstype,M0,pBasis)
% estimate a PD from the data
%  Guesstype=1 or T1reg=0 use the sum of sqr (or M0) as PD
%  Guesstype=2 randum start
% Guesstype=3 start or T1reg>0  with the PD R1 relation as describe in the literature.
% Guesstype=4 use ridge regration Bilinear fit to start.
% AM  & BW VISTASOFT Team, 2013

if notDefined('Guesstype')
Guesstype=0;
end

if (T1reg==0 || Guesstype==1) % sum of squr
PDinit = sqrt(sum(M0.^2,2));  
end

if  Guesstype==2 %rand PD
PDinit = rand(size(M0,1),1);
end

if (T1reg>0 || Guesstype==3)% T1 PD linear relations
    
          PDinit=1./(R1*0.42+0.95); %this is the teortical T1 PD relationship see reference at Mezer et. al 2013
end

if (Guesstype==4)
    %ridge regration
   maxLoops=200;
 sCriterion = 1e-3;  % Stopping criterion
 lambda=1;
 BLFit_RidgeReg = pdBiLinearFit(M0, pBasis, ...
                 lambda, maxLoops, sCriterion, [], 0 );
PDinit=BLFit_RidgeReg.PD;
end



%guess the g0

if ~notDefined('pBasis')
 
    nVoxels=size(M0,1);
    Ncoils=size(M0,2);
    nPolyCoef=size(pBasis,2);
    
     G  = zeros(nVoxels,Ncoils);
    g0 = zeros(nPolyCoef,Ncoils);
    for ii=1:Ncoils
        G(:,ii)  = M0(:,ii) ./ PDinit(:);         % Raw estimate
        g0(:,ii) = pBasis(:,:) \ G(:,ii);  % Polynomial approximation
    end
else
    g0=[];
end