function err = errFitNestBiLinearTissueT1reg_full(g,M0,pBasis,nPositions,nCoils,R1basis,RegWeight,TissueMask,M0mask,T)
% Bilinear estimation subject to tisuue spesipic T1 regularization as well
%
%  err =
%   errFitNestBiLinearTissueT1reg_full(g,M0,pBasis,nPositions,nCoils,R1basis,RegWeight,TissueMask,M0mask)
%
% AM/BW(c) VISTASOFT Team, 2013

%estimate coil coefficients across the volume
G = pBasis*g;

% Estimate the best PD for each position a linear sulotion
% This makes it a nested bilinear problem

PD = zeros(nPositions,1);
for ii=1:nPositions
    %use=M0mask(ii,:);
    PD(ii) = G(ii,M0mask(ii,:))' \ M0(ii,M0mask(ii,:))';
end


% Given the known T1, which is sent in, we have an expectation that it will
% be linearly related to the estimated PD via this equation:
%
%  1/PD = C1/T1+ C2
%
% So, we have a PD predicted from the T1.  These are the coefficients of
% the estimated linear relationship
% co     = R1basis \ ( 1./PD(:) );
% PDpred = R1basis*co;
%T= (unique(TissueMask));
PDpred=zeros(size(PD));
for ii=1:length(T); 
    if ii>0
        mask=find(TissueMask==T(ii));
        PDpred(mask)     = R1basis(mask,:)* (R1basis(mask,:) \ ( 1./PD(mask) ));
    end
end
PDpred     = 1 ./ PDpred;

% Normalize by the first PD value
% G  = G  .* PD(1);
% PD = PD ./ PD(1);

G  = G  .* mean(PD(TissueMask>0));
PDpred = PDpred ./ mean(PD(TissueMask>0));
PD = PD ./ mean(PD(TissueMask>0));

% get the predicted M0 for all of the coils
M0P = G.*repmat( PD,1,nCoils);


% mrvNewGraphWin; plot(PD(:),PDpred(:),'o')
% PDpred     = R1basis* (R1basis \ ( 1./PD(:) ));
% PDpred     = 1 ./ PDpred;

% The error is PD - PDpred
%
% An alternative is to calculate the correlation coefficient
% c = corrcoef(R1basis(:,2), 1./PD(:)); err = RegWeight*(1 - c(1,2));
% And then use this number at the end

% The error is a vector with positive and negative values representing the
% M0 difference and the T1 linearity failures
%err = [ M0(M0mask) - M0P(M0mask); (nCoils*RegWeight)*(PD(TissueMask>0) - PDpred(TissueMask>0))];
err = [ M0(M0mask) - M0P(M0mask); (RegWeight)*(PD(TissueMask>0) - PDpred(TissueMask>0))];

end