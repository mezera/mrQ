function [params, STR] = polyGetPhantomCoef(Ncoils,dim,sample,whichCoils)
% Get polynomial coefficients from the phantom data
%
%  [params STR] = polyGetPhantomCoef(Ncoils,dim,sample,whichCoils)
%
% Load the polynomial coefficents defining the coil gains estimated from a
% phantom on 32 channel Nova coil at three different spatial locations
%
% Ncoils: Number of coil parameters. if diferent the whichCoils it will be
%         the number of whichCoils.  (default: 2)
% dim:          1D 2D 3D poly coefisents  (default: 2D)
% sample:       what sample from the phantom to use (1,2,3). (default: 1)
% whichCoils:   which coil sample to use. a vector with number between 1 to
%               32. (default: 1:Ncoils) 
% Output:   
%   params: The polynomial coefficients
%   STR:   Strings defining the polynomial terms (X^2) and so forth
%
% Example:
%   Ncoils = 32; dim = 2; sample = 1;
%   [params STR] = polyGetPhantomCoef(Ncoils,dim,sample)
%
% AM Copyright VISTASOFT Team, 2013

% Check input parameters
if notDefined('Ncoils'), Ncoils=2; end
if notDefined('sample'), sample=1; end
if notDefined('whichCoils'), whichCoils=1:Ncoils; end
if notDefined('dim'), dim=2; end

if dim == 1,      ed = 3;
elseif dim == 2,  ed =6;
elseif dim == 3,  ed=10;
end

% Check the number of coils
if Ncoils ~= length(whichCoils)
    Ncoils = length(whichCoils);
end

% Load coil gains into the parameter C
load('CoilGains')

% Get the values for a spatial sample, some number of polynomial dimensions
% and some list of coils
params(1:ed,1:Ncoils) = C{sample}.params(1:ed,whichCoils);

% Not sure what this is.
STR = C{sample}.STR(1:ed);

return
