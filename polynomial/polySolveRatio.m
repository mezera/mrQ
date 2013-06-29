function eGains = polySolveRatio(polyRatioMat,M0pairs)
% Solve for the gain parameters of the polynomial coefficient
%
%  eGains = polySolveRatio(polyRatioMat,M0pairss)
%
% Inputs:
%  polyRatioMat:  The large matrix such that 0 = polyRatioMat * gains
%                 See polyCreateRatio
%  M0pair:        The pair of measurements corresponding to each row
%
% Returns
%  eGains:  Estimated coil gain coefficients
%
% BW, AM Copyright Vistasoft Team, 2013


%% The rows of polyRatioMat are ratios of M0 measurements
%
%  These measurements are from two coils at a voxel.
%  We could scale the row to reflect the SNR of that ratio.
%  For example, if one of the M0 values is small, that would reduce our
%  confidence in the ratio.
%
% Note, we have a lot of coil and voxels.  So we can get rid of quite a few
% (even zeroing them out) and still solve this equation.

% if exist('M0pairs','var')
%     wgts = pairs2weights(M0pairs);
%     polyRatioMat = diag(wgts)*polyRatioMat;
% end

%% I am concerned that the gain coefficients are so different in size
%
% The DC coefficient is 1000 times bigger than the other ones.  
% In fact, they decline precipitously across order.  So, maybe we should be
% doing something to balance these better.
%% We could re-write the M0 matrix to represent virtual coils
%
%  Each column of M0 is a real coil.  If we analyze the 

% Solve for the eigenvectors of polyRatioMat.  The first one has the
% smallest singular value and is thus the best estimate of the solution to
% the equation 0 = polyRatioMat * g
%
%  
[U, ~] = eig(polyRatioMat'*polyRatioMat);
 
% % We scale the gain vector so that the first entry of the first coil
% % estimate (the constant) is 1.
eGains = U(:,1)/U(1,1);


end

%%  Comments and notes - one way to do it
%
% This is an equivalent way to calculate this
% col = polyRatioMat(:,1);
% col = -1*col; 
% 
% reducedMat = polyRatioMat(:,2:end);
% nCols = size(reducedMat,2);
% 
% % One possible solution
% %  col = reducedMat*g;
% g = reducedMat\col;
% %  reducedMat*g - Should equal
% eGains2 = [1 ; g];
% % polyRatioMat*eGains2 should be 0
% % polyRatioMat*eGains2;
% 
% % mrvNewGraphWin; plot(eGains(:),eGains2(:),'.')

%% More comments
% % This is trying to control the singular values
% %
% %
% %
% 
% [U, S, V] = svd(reducedMat);
% s = diag(S);
% s = 1 ./ s;
% zeroSingularValues = 26;
% s(zeroSingularValues:end) = 0;
% Sr = diag(s);
% U = U'; U = U(1:nCols,:);
% p = V * Sr * U;  % This is the pseudo inverse
% 
% % Now solve
% eGains3 = [1; p*col];
% % mrvNewGraphWin; plot(eGains(:),eGains3(:),'.')
%
%
%%


function wgts = pairs2weights(M0pairs)

%% Default and does nothing
wgts = ones(size(M0pairs,1),1);

%% An experiment
% min(M0pairs,[],2);
% wgts = min(M0pairs,[],2);
% wgts = wgts.^2;

%% ANother experiment
% top1 = prctile(M0pairs(:,1),50);
% top2 = prctile(M0pairs(:,2),50);
% lst1 = (M0pairs(:,1) > top1);
% lst2 = (M0pairs(:,2) > top2);
% lst = lst1 .* lst2;
% wgts = zeros(size(M0pairs,1),1);
% wgts(logical(lst)) = 1;
% 
% sum(wgts(:))

end


