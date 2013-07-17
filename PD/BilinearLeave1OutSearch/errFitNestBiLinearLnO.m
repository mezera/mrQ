function err=errFitNestBiLinearLnO(g,M0,pBasis,nPositions,Fcoils,Tcoils)
% err=errFitNestBiLinearLnO(g,M0,pBasis,nPositions,Fcoils,Tcoils)
% we fit PD with some coils (Fcoils) and calculate the error of Mo prediction on other coils
% (Tcoils). this make a leave N out fitting. that might be less saseptibe
% to noise.
err=[];
%estimate coil coefisets
G = pBasis*g;

% get the best PD for each position a linear sulotion
%this mske  it a nested biliner problem
for jj=1:length(Fcoils)
PD = zeros(nPositions,1);
    nCoils=length(Tcoils{jj});
    for ii=1:nPositions
        PD(ii) = G(ii,Fcoils{jj})' \ M0(ii,Fcoils{jj})';
    end
    %normelize by the first PD value
    G  = G  .* PD(1);
    PD = PD ./ PD(1);
    
    % get the predigted M0 from the Test coils
    M0P = G(:,Tcoils{jj}).*repmat( PD,1,nCoils);
    
    
    errT=M0(:,Tcoils{jj})-M0P;
    err=[err(:) ;errT(:)];
end
end

