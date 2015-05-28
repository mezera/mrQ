function [T1Est, bMagEst, aMagEst, res] = ...
  T1FitNLSPR(data, nlsS)
    
% [T1Est, bMagEst, aMagEst, res] = ...
%  T1FitNLSPR(data, nlsS)
%
% Finds estimates of T1, a, and b using a nonlinear least
% squares approach together with phase restoration. 
% The model +-|aMag + bMag*exp(-t/T1)| is used. 
% The residual is the rms error between the data and the fit. 
% 
% INPUT:
% data - the absolute data to estimate from
% nlsS - struct containing the NLS search parameters and
%        the data model to use
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 


if ( length(data) ~= nlsS.N )
  error('nlsS.N and data must be of equal length!')
end

%Make sure the data comes in increasing TI-order
[tVec,order] = sort(nlsS.tVec); 
data = squeeze(data); 
data = data(order);
theExp = nlsS.theExp(order,:);

aEstTmp = zeros(1,nlsS.N);
bEstTmp = zeros(1,nlsS.N);
T1EstTmp = zeros(1,nlsS.N);
resTmp = zeros(1,nlsS.N);

% Make sure data is a column vector
data = data(:);

switch(nlsS.nlsAlg)
  case 'grid'
    try
      nbrOfZoom = nlsS.nbrOfZoom;
    catch
      nbrOfZoom = 0; % No zoom
    end

    for ii = 1:nlsS.N
      %Setting (ii-1)th first elements to minus
      dataTmp = data.*[-ones(ii-1,1); ones(nlsS.N - (ii-1),1)];

      % The sum of the data
      ySum = sum(dataTmp);

      switch(nlsS.dataModel)
        case{'ab'}
          rhoTyVec = (dataTmp.'*theExp).' - ...
            1/nlsS.N*sum(theExp,1)'*ySum;
      end

      rhoNormVec = nlsS.rhoNormVec;
      % The maximizing criterion
      [tmp,ind] = max( abs(rhoTyVec).^2./rhoNormVec );

      T1Vec = nlsS.T1Vec; % Initialize the variable
      if nbrOfZoom > 1 % Do zoomed search
        try
          T1LenZ = nlsS.T1LenZ; % For the zoomed search
        catch
          T1LenZ = 21; % For the zoomed search
        end
        for k = 2:nbrOfZoom
          if( ind > 1 && ind < length(T1Vec) )
            T1Vec = linspace(T1Vec(ind-1),T1Vec(ind+1),T1LenZ)';
          elseif(ind == 1)
            T1Vec = linspace(T1Vec(ind),T1Vec(ind+2),T1LenZ)';
          else
            T1Vec = linspace(T1Vec(ind-2),T1Vec(ind),T1LenZ)';
          end
          alphaVec = 1./T1Vec;
          theExp = exp( -tVec*alphaVec' );
          yExpSum = (dataTmp.'*theExp).';

          switch(nlsS.dataModel)
            case{'ab'}
              rhoNormVec = ...
                sum( theExp.^2, 1)' - ...
                1/nlsS.N*(sum(theExp,1)').^2;
              rhoTyVec = yExpSum - ...
                1/nlsS.N*sum(theExp,1)'*ySum;
          end

          [tmp,ind] = max( abs(rhoTyVec).^2./rhoNormVec );
        end
      end %of zoom

      % The estimated parameters
      T1EstTmp(ii) = T1Vec(ind);
      bEstTmp(ii) = rhoTyVec(ind)/ rhoNormVec(ind);
      switch(nlsS.dataModel)
        case{'ab'}
          aEstTmp(ii) = 1/nlsS.N*(ySum - bEstTmp(ii)*sum(theExp(:,ind)));
      end
      
      % Compute the residual
      switch(nlsS.dataModel)
        case{'ab'}
          modelValue = aEstTmp(ii) + bEstTmp(ii)*exp(-tVec/T1EstTmp(ii));
      end

      resTmp(ii) = 1/sqrt(nlsS.N) * norm(1 - modelValue./dataTmp);
    end %of for loop

otherwise
  error('Unknown search method!')
end

[res,ind] = min(resTmp);
aMagEst = aEstTmp(ind);
bMagEst = bEstTmp(ind);
T1Est = T1EstTmp(ind);
