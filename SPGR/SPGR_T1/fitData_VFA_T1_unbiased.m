function [t1Unbiased t1Biased] = fitData_VFA_T1_unbiased (data, flipAngles, TR, b1Map)

% function [t1Unbiased t1Biased] = fitData_VFA_T1_unbiased (t1Data, flipAngles, b1Map)
% -----------------------------------------------------------
% This function performs a weighted-least squares data fit
% on SPGR T1 data set
% INPUTS:
% t1Data: width x length x slices x flipAngles matrix
% flipAngles: a vector of flip angles (in degrees) that corresponds to the
% t1Data matrix size
% TR: in ms
% b1Map: a width x length x slices matrix that contains the relative flip
% angle (i.e. if nominal alpha is 60, and measured alpha is 61, then b1Map
% = 61/60
%
% Nikola Stikov 2013
%
if (nargin < 4)
    dataSize = size(data);
    b1Map = ones(dataSize(1:end-1));
end

dims = size(data);
t1Biased = zeros(dims(1:end-1));
t1Unbiased = zeros(dims(1:end-1));

for ii=1:dims(1)
    ii;
    for jj=1:dims(2)
        for kk=1:dims(3)
            if (sum(data(ii, jj, kk, :))==0 ||  isnan(sum(data(ii, jj, kk, :))) ||   isinf(sum(data(ii, jj, kk, :))));
                
            else
                y = squeeze(data(ii, jj, kk, :))./sin(flipAngles/180*pi*b1Map(ii, jj, kk))';
                x = squeeze(data(ii, jj, kk, :))./tan(flipAngles/180*pi*b1Map(ii, jj, kk))';
                    test = polyfit(x, y, 1);
                    slopeBiased = test(1);
                    result = abs(-TR./log(slopeBiased));
                    if (~isnan(result))% && result < 5)
                        t1Biased(ii,jj,kk) = result;
                    end
                    
                    % show we nclude B1 in the flip angle for theweights as well??
                    weights = (sin(flipAngles/180*pi)./(1 - slopeBiased.*cos(flipAngles/180*pi))).^2;
                    weights(isinf(weights)) = 0; %remove points with infinite weight
                    
                    test2 = polyfitweighted(x, y, 1, weights');
                    slopeUnbiased = test2(1);
                    result = abs(-TR./log(slopeUnbiased));
                    if (~isnan(result) )%&& result < 5)
                        t1Unbiased(ii,jj,kk) = result;
                    end
                        
                end
            end
        end
    end