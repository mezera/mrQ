function [params STR]=GetPantomPolyCoef(Ncoils,dim,sample,whichCoils)
%[params STR]=GetPantomPolyCoef(Ncoils,dim,sample,whichCoils)
%
%the function load the poly coefficent as estimate on a phantom on 32
%chancel nova coil in three differnent locations in space
%Input:
% Ncoils           number of coil parameters. if diferent the whichCoils it will be
%                      the number of whichCoils.  (defult 2)
% dim              1D 2D 3D poly coefisents. defult 2D
%sample          what sample from the phantom to use (1,2,3). defult 1
%whichCoils   which coil sample to use. a vector with number between 1 to 32. (defult:   1:Ncoils)
%
% Output:   params the poly coefisents that are relate to STR polyinomyal

% define undefine input to defult
if notDefined('Ncoils')
    Ncoils=2;
end

if notDefined('sample')
    sample=1;
end
if notDefined('whichCoils')
    whichCoils=1:Ncoils;
end

if notDefined('dim')
    dim=2;
end
if dim==1
    ed=3;
elseif dim ==2
    ed =6;
elseif dim==3
    ed=10;
end

%check the number of coils
if Ncoils~=length(whichCoils)
    Ncoils=length(whichCoils);
end

% load data
load ('CoilGains')

%get the values
params(1:ed,1:Ncoils)=C{sample}.params(1:ed,whichCoils);
STR=C{sample}.STR(1:ed);

