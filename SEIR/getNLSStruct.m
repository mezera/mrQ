function nlsS = getNLSStruct_v3( extra, dispOn, zoom)

% extra.tVec    : defining TIs 
%                 (not called TIVec since it looks too much like T1Vec)
% extra.T1Vec   : defining T1s
% dispOn        : 1 - display the struct at the end
%                 0 (or omitted) - no display
% zoom          : 1 - do a non-zoomed search
%                 2 - do a zoomed search (slower search)
% Data Model    : a + b*exp(-TI/T1) 
%    
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University    

%%

nlsS.tVec    = extra.tVec(:);
nlsS.N       = length(nlsS.tVec);
nlsS.T1Vec   = extra.T1Vec(:);
nlsS.T1Start = nlsS.T1Vec(1);
nlsS.T1Stop  = nlsS.T1Vec(end);
nlsS.T1Len   = length(nlsS.T1Vec);

% The search algorithm to be used
nlsS.nlsAlg = 'grid'; % Grid search

% Display the struct so that the user can see it went ok
if nargin < 2 
  dispOn = 0;
end

% Set the number of times you zoom the grid search in, 1 = no zoom
% Setting this greater than 1 will reduce the step size in the grid search
% (and slow down the fit significantly)

if nargin < 3
	nlsS.nbrOfZoom = 1;
else
	nlsS.nbrOfZoom = zoom;
end

if nlsS.nbrOfZoom > 1
    nlsS.T1LenZ = 21; % Length of the zoomed search
end
		
nlsS.dataModel = 'ab';

if dispOn
  nlsS
end

% Set the help variables that can be precomputed
switch(nlsS.nlsAlg)
  case{'grid'}
    alphaVec = 1./nlsS.T1Vec; 
    nlsS.theExp = exp( -nlsS.tVec*alphaVec' );
    nlsS.rhoNormVec = ...
        sum( nlsS.theExp.^2, 1)' - ...
        1/nlsS.N*(sum(nlsS.theExp,1)').^2;    
end 

