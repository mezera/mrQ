function [T1Hat, bHat, aHat, residual] = T1FitNLS(data, nlsS)
    
% [T1Hat, bHat, aHat, residual] = T1FitNLS(data, nlsS)
%
% Finds estimates of T1, a, and b using a nonlinear least
% squares approach. The model a+b*exp(-t/T1) is used. 
% The residual is the rms error between the data and the fit. 
% 
% INPUT:
% data - the data to estimate from
% nlsS - struct containing the NLS search parameters and
%        the data model to use
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 


if ( length(data) ~= nlsS.N )
  error('nlsS.N and data must be of equal length!')
end

aHat = 0;

% Make sure data is a column vector
data = data(:);

switch(nlsS.nlsAlg)
  case 'grid'
    try
      nbrOfZoom = nlsS.nbrOfZoom;
    catch
      nbrOfZoom = 0; % No zoom
    end
    
    % The sum of the data
    ySum = sum(data);
    
    switch(nlsS.dataModel)
      case{'ab'}
        rhoTyVec = (data.'*nlsS.theExp).' - ...
          1/nlsS.N*sum(nlsS.theExp,1)'*ySum;
        theExp = nlsS.theExp;
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
        theExp = exp( -nlsS.tVec*alphaVec' );
        yExpSum = (data.'*theExp).';
        
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
    end
    
    % The estimated parameters
    T1Hat = T1Vec(ind);
    bHat = rhoTyVec(ind)/ rhoNormVec(ind);
    switch(nlsS.dataModel)
      case{'ab'}
        aHat = 1/nlsS.N*(ySum - bHat*sum(theExp(:,ind)));
    end
    
  otherwise
    error('Unknown search method!')
end

% Compute the residual
switch(nlsS.dataModel)
  case{'ab'}
    modelValue = aHat + bHat*exp(-nlsS.tVec/T1Hat);
end

residual = 1/sqrt(nlsS.N) * norm(1 - modelValue./data);
