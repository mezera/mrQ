function [transmit, receive_A, receive_D] = getGainFromDcm(dcmPath)
% Get transmit and receive gain from a dicom header
%
% [transmit, receive_A, receive_D] = getGainFromDcm(dcmPath)
%
% Inputs:
% dcmPath   - Path to a dicom image
%
% Outputs:
% transmit  - Transmit gain (.1 db)
% receive_A - Receive gain analog (.1db)a
% receive_D - Receive gain digital (.1 db)
%
% Written by Jason D. Yeatman April 4 2012

if notDefined('dcmPath') || ~exist(dcmPath,'file')
    error('\nPlease give a valid path to a dicom image')
end

d = dicominfo(dcmPath);

transmit  = d.Private_0019_10f9;
receive_A = d.Private_0019_108a;
receive_D = d.Private_0019_108b;