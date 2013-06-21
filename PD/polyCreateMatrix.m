function [pMatrix,s] = polyCreateMatrix(nSamples,order,dimension)
% Build 2D polynomial matrix.  Will generalize to other cases later
%
%    [pMatrix,s] = polyCreateMatrix(nSamples,order,dimension)
%
% pMatrix:  Polynomial basis functions for nth order and some number of
% dimensions.   This matrix does NOT include the constant
% s:     Spatial samples
%
% Example:
%   nSamples = 20; order = 2; dimension = 2;
%   pMatrix = polyCreateMatrix(20,2,2);
%   imagesc(pMatrix);
%
%  [M2, O] = polyCreateMatrix(20,2,2);
%  M2 = [O,M2];
%  imagesc(pMatrix);
%
% BW Copyright vistasoft 2013

s = -nSamples:nSamples;

switch order
    case 2  % 2nd order polynomial
        if dimension == 1
            X = s(:); 
            X2 = X(:).^2;
            pMatrix = [ones(size(X)), X, X2];
        elseif dimension == 2
            [X, Y] = meshgrid(s,s);
            X = X(:);
            Y = Y(:);
            X2 = X(:).^2;
            Y2 = Y(:).^2;
            XY = X(:).*Y(:);
            % X, X^2, Y, Y^2, XY
            pMatrix = [ones(size(X)), X, X2, Y, Y2, XY];
        else
            error('Not yet implemented')
        end
        
    otherwise
        error('Order %d not built',order);
end

% Sometimes the user will want a column of ones of the proper size.
if nargout > 1
    O = ones(length(s)^2,1);
end

return
