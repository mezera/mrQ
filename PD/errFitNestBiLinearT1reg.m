function err=errFitNestBiLinearT1reg(g,M0,pBasis,nPositions,nCoils,R1basiss,RegWight)


%estimate coil coefisets
G = pBasis*g;

% get the best PD for each position a linear sulotion
%this mske  it a nested biliner problem
PD = zeros(nPositions,1);
for ii=1:nPositions
    PD(ii) = G(ii,:)' \ M0(ii,:)';
end
%normelize by the first PD value
G  = G  .* PD(1);
PD = PD ./ PD(1);

% get the predigted M0
M0P = G.*repmat( PD,1,nCoils);

% solve ls problem for linear relations  1/PD= C1/T1+ C2
co=R1basiss\(1./PD(:));
PDpred=R1basiss*co;
err=[M0-M0P; RegWight(PD-PDpred)];
end