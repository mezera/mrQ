function [pMatrix,s] = polyCreateMatrix(nSamples,order,dimension)
% Build 2D polynomial matrix 
%
%    [pMatrix,s] = polyCreateMatrix(nSamples,order,dimension)
%
% nSamples: Runs from -nSamples to +nSamples
% order   : polynomial order (linear, quadratic, cubic)
% dimension: 1, 2, or 3
%
% pMatrix: Polynomial basis functions for nth order and some number of
%          dimensions.   This matrix does NOT include the constant
% s:       Spatial samples
%
% Currently implemented for 2nd order, 1D, 2D and 3D
%
% Example:
%   nSamples = 20; order = 2; dimension = 2;
%   pMatrix = polyCreateMatrix(20,2,2);
%   imagesc(pMatrix);
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
            X = X(:);     Y = Y(:);
            X2 = X(:).^2; Y2 = Y(:).^2;
            XY = X(:).*Y(:);
            pMatrix = [ones(size(X)), X, X2, Y, Y2, XY];
            
         elseif dimension == 3
            [X, Y, Z] = meshgrid(s,s,s);
            X = X(:); Y = Y(:); Z = Z(:);
            X2 = X(:).^2; Y2 = Y(:).^2; Z2 = Z(:).^2;
            XY = X(:).*Y(:); XZ = X(:).*Z(:);  YZ = Y(:).*Z(:);
            pMatrix = [ones(size(X)), X, X2, Y, Y2, XY  Z  Z2 XZ YZ];
            
        else
            error('Not yet implemented')
        end
        
    otherwise
        error('Order %d not built',order);
end

return
