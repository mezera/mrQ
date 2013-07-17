function err=errFitNestBiLinear(g,M0,pBasis,nPositions,nCoils)

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


err=M0-M0P;
end