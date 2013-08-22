function rootPath=mrqRootPath()
% Return root of the mrQ directory
% 
%        rootPath = mrqRootPath;
%
% This function MUST reside in the directory at the base of the mrQ
% directory structure
%
% Wandell Copyright Vistasoft Team, 2013

rootPath = fileparts(which(mfilename));

return
