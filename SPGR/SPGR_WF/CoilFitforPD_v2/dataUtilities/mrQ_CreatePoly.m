function [pBasis, pTerms, nSamples] = ...
    mrQ_CreatePoly(boxS,pOrder,sDim,BasisFlag)
% Build 2D polynomial matrix
%
%     [pBasis, pTerms, nSamples] = ...
 %        mrQ_CreatePoly(boxS,pOrder,sDim,BasisFlag)
%
% boxS: the     X Y Z size of the box
% pOrder  :     polynomial order (linear, quadratic, cubic)
% sDim:    :    1, 2, or 3
% BasisFlag:  Normalize pBasis if true to be unit lenght. if 'svd' use svd
%                   to Orthogonalize pBasis. if 'qr' use qr to  Orthogonalize thepBasis.
%
% pBasis: Polynomial basis functions for nth order and some number of
%          dimensions.   This matrix does NOT include the constant
% spatialSamples:       Spatial samples
% pTerms:  String defining the polynomial terms
%
% Currently implemented for 2nd order, 1D, 2D and 3D
%
% Example:
%   boxS = [12 10]; order = 3; dimension = 2;
%   pBasis = polyCreateMatrix(boxS,order,dimension);
%   imagesc(pBasis);
%
% AM & BW Copyright vistasoft 2013

%%
if notDefined('pOrder'), pOrder = 3; end
if notDefined('sDim'), sDim = 3; end
if notDefined('BasisFlag'), BasisFlag = false; end

nSamples=floor(boxS/2);
spatialSamplesX = -nSamples(1):nSamples(1);
spatialSamplesY = -nSamples(2):nSamples(2);
spatialSamplesZ = -nSamples(3):nSamples(3);

%%
switch pOrder
    case 1  % 1st order polynomial
        if sDim == 1
            X = spatialSamplesX(:);
            pBasis = [ones(size(X)), X];
            pTerms = '[1, X]';
         %   W=[0 1];
        elseif sDim == 2
            [X, Y] = meshgrid(spatialSamplesX,spatialSamplesY);
            X = X(:);     Y = Y(:);
            pBasis = [ones(size(X)), X, Y];
            pTerms = '[1, X, Y]';
         %   W=[0 1 1];
        elseif sDim == 3
            [X, Y, Z] = meshgrid(spatialSamplesX,spatialSamplesY,spatialSamplesZ);
            X = X(:); Y = Y(:); Z = Z(:);
            pBasis = [ones(size(X)), X, Y  Z ];
            pTerms = '[1, X, Y  Z ]';
            % W=[0 1 1 1];
        else
            error('Not yet implemented')
        end
        
    case 2  % 2nd order polynomial
        if sDim == 1
            X = spatialSamplesX(:);
            X2 = X(:).^2;
            pBasis = [ones(size(X)), X, X2];
            pTerms  = '[1, X, X2]';
         %   W=[0 0.1 1];
            
        elseif sDim == 2
            [X, Y] = meshgrid(spatialSamplesX,spatialSamplesY);
            X = X(:);     Y = Y(:);
            X2 = X(:).^2; Y2 = Y(:).^2;
            XY = X(:).*Y(:);
            % Six parameters
            pBasis = [ones(size(X)), X, X2, Y, Y2, XY];
            pTerms  = '[1, X, X2, Y, Y2, XY]';
      %      W=[0 0.1 1 0.1 1 1];
        elseif sDim == 3
            [X, Y, Z] = meshgrid(spatialSamplesX,spatialSamplesY,spatialSamplesZ);
            X = X(:); Y = Y(:); Z = Z(:);
            X2 = X(:).^2; Y2 = Y(:).^2; Z2 = Z(:).^2;
            XY = X(:).*Y(:); XZ = X(:).*Z(:);  YZ = Y(:).*Z(:);
            % Ten parameters
            pBasis = [ones(size(X)), X, X2, Y, Y2, XY  Z  Z2 XZ YZ];
            pTerms  = '[1, X, X2, Y, Y2, XY  Z  Z2 XZ YZ]';
   %           W=           [0 0.1 1   0.1  1   1    0.1 1 1 1];
        else
            error('Not yet implemented')
        end
        
    case 3
        
        if sDim == 1
            X = spatialSamplesX(:);
            X2 = X(:).^2;
            X3 = X(:).^3;
            pBasis = [ones(size(X)), X, X2, X3];
            pTerms  = '[1, X, X2, X3]';
  %               W=           [0  0.01  0.1       1    ];
        elseif sDim == 2
            [X, Y] = meshgrid(spatialSamplesX,spatialSamplesY);
            X = X(:);     Y = Y(:);
            X2 = X(:).^2; Y2 = Y(:).^2;
            X3 = X(:).^3; Y3 = Y(:).^3;
            XY = X(:).*Y(:);
            X2Y = X2.*Y(:); XY2 = X(:).*Y2(:);
            
            % 10 parameters
            pBasis = [ones(size(X)), X, X2, Y, Y2, XY, X3, Y3, X2Y, XY2];
            pTerms  = '[1,  X,     X2,   Y,    Y2,  XY, X3, Y3, X2Y, XY2]';
%             W=            [0 0.01  0.1  0.01  0.1  0.1    1      1     1        1     ];

        elseif sDim == 3
            [X, Y, Z] = meshgrid(spatialSamplesX,spatialSamplesY,spatialSamplesZ);
            X = X(:); Y = Y(:); Z = Z(:);
            X2 = X(:).^2; Y2 = Y(:).^2; Z2 = Z(:).^2;
            X3 = X(:).^3; Y3 = Y(:).^3; Z3 = Z(:).^3;
            
            XY = X(:).*Y(:); XZ = X(:).*Z(:);  YZ = Y(:).*Z(:);
            XYZ = X(:).*Y(:).*Z(:);
            X2Y = X2.*Y(:); X2Z = X2.*Z(:);
            XY2 = X(:).*Y2(:); XZ2 = X(:).*Z2(:);
            Y2Z = Y2(:).*Z(:); YZ2 = Y(:).*Z2(:);
            
            % Twenty parameters
            pBasis = [ones(size(X)), X, X2, Y, Y2, XY ,X3, Y3, X2Y, XY2, Z, Z2, Z3, XZ, YZ, XYZ, X2Z, Y2Z, XZ2, YZ2];
            pTerms  = '[1, X,    X2,   Y,   Y2, XY ,X3, Y3, X2Y, XY2, Z,      Z2,    Z3,  XZ, YZ, XYZ, X2Z, Y2Z, XZ2, YZ2]';
   %         W=             [0 0.01 0.1 0.01 0.1  0.1   1      1     1          1    0.01  0.1    1      0.1  0.1  1         1        1        1       1 ] ;
        end
        
    otherwise
        error('Order %d not built',pOrder);
end

if BasisFlag
    if ischar(BasisFlag)
        BasisFlag = mrvParamFormat(BasisFlag);
        switch BasisFlag
            case {'svd'}
                % SVD Orthogonalize the pBasis
                nCols = size(pBasis,2);
                [U, ~, ~] = svd(pBasis);
                pBasis = U(:,1:nCols);
            case {'qr'}
                % QR Orthogonalize the pBasis
                nCols = size(pBasis,2);
                [Q R] = qr(pBasis);
                pBasis = Q(:,1:nCols);
                % Make sure the first one (mean) is positive
                if pBasis(1,1) < 0, pBasis(:,1) = -1*pBasis(:,1); end
            otherwise
                error('BasseFlag %d not built',BasseFlag);
        end
    else
        % Adjust the length of basis vectors, but don't orthogonalize.
        sFactor = sqrt(diag(pBasis'*pBasis));
        pBasis = pBasis * diag(1./sFactor);
    end
end

return
