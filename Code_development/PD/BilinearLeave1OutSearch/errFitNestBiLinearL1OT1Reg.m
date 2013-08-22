function err=errFitNestBiLinearL1OT1Reg(g,M0,pBasis,nPositions,Fcoils,Tcoils,R1basis,RegWeight)
%err=errFitNestBiLinearL1OT1Reg(g,M0,pBasis,nPositions,Fcoils,Tcoils,,R1basis,RegWeight)
% we fit PD with some coils (Fcoils) and calculate the error of Mo prediction on other coils
% (Tcoils). this make a leave N out fitting. that might be less saseptibe
% to noise.
%last the Bilinear estimation subject to T1 regularization as well (linear
%relation to T1)

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

% now get the PD according to all coils 

PD = zeros(nPositions,1);
for ii=1:nPositions
    PD(ii) = G(ii,:)' \ M0(ii,:)';
end

PDpred     = R1basis* (R1basis \ ( 1./PD(:) ));
PDpred     = 1 ./ PDpred;

%join the errors
err = [ err(:); RegWeight*(PD(:) - PDpred(:))];

end

