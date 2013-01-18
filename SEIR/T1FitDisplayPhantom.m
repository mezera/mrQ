function roi = T1FitDisplayPhantom(...
  loadStr, saveStr, method, roi)

% T1FitDisplayPhantom(loadStr, saveStr, method)
% or
% T1FitDisplayPhantom(loadStr, saveStr, method, roi)
%  
% Mandatory:
%    loadStr: the data that the fit was performed on, to be loaded
%    saveStr: the results of the fit to be loaded
%    method:  what fitting method to use,
%             NLS, LM, LMsph, LMsphMag, NLS_mag or LM_mag
%
% Optional (if omitted, the ROI is manually chosen): 
%    roi - preset of ROI
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University

% Load the original data 
load(loadStr);

dims = size(data);

if numel(dims) > 3
  nbslice = dims(3); % See how many slices there are
else
  nbslice = 1;
  tmpData(:,:,1,:) = data; % Make data a 4-D array regardless of number of slices
  data = tmpData;
  clear tmpData;
end

load(saveStr); % This loads ll_T1
T1_original = ll_T1; % Store it so ll_T1 can be a single slice;
Mask_original = mask; % Store it so mask can be a single slice

% Inserting a short pause, otherwise some computers seem
% to get problems
pause(1);

if nargin == 3
  nbrois = input('How many ROIs?   --- ');
else
  nbrois = size(roi,3);
  roiOrg = roi;
end

% Display T1 data
showMontage(squeeze(abs(ll_T1(:,:,:,1))),[],'colorMap');
axis equal
axis off
colorbar
title('T1 values for all slices');
customFormat

for r = 1:nbrois
  while(true)
    if (nbslice ~= 1)
      zz = input(['Slice number for ROI ' num2str(r)...
        '? (Enter number 1 to ' num2str(nbslice) ')   --- '], 's');
      zz = cast(str2num(zz), 'int16');
      if (isinteger(zz) & zz >0 & zz <= nbslice) break;
      end
    else
      zz = 1; break;
    end
  end
  datali = squeeze(data(:,:,zz,end));
  ll_T1 = squeeze(T1_original(:,:,zz,:));
  T1 = ll_T1(:,:,1);
  mask = squeeze(Mask_original(:,:,zz));
  figure(111),
  imshow(abs(datali),[])
  
  if (nargin==3)
	  figure(111)
	  % Select ROI
	  roiname = input(['ROI ' num2str(r) ' name?   --- '], 's');
	  disp('Draw the ROI, double click inside when done')
	  roi = roipoly;
	  close(111)
  else
    roiname = 'ROI';
    roi = roiOrg(:,:,r);
  end

  % Mask out points where grid search hit boundary
  mask(find(T1==extra.T1Vec(1) | T1==extra.T1Vec(end))) = 0;
  T1 = T1.*mask;
  T1roi = T1.*roi;
  % Check if all elements are zero 
  if all( all(T1roi == 0) )
    error('All T1 estimates in the ROI is zero, no histogram can be plotted')
  end
    
  x = min(T1roi(find(T1roi>0))):1:max(T1roi(find(T1roi>0)));
  x = x';
  h = histc(T1roi(find(T1roi>0)),x); % centered on integers
  hs = h;  
  [sigma,mu,A] = customGaussFit(x,hs);
  figure
  bar(x,hs,1) % centered on integers
  hold
  plot(x, A*exp(-(x-mu).^2/(2*sigma^2)),'.r')
  
  titlename = [roiname, ': ', num2str(sum(h)), ' pixels, mode ', num2str(round(mu)),...
	  ' ms, \sigma = ', num2str(round(sigma))]
  title(titlename)
  ylabel('Number of pixels')
  xlabel('T_1 [ms]')
  customFormat
  clear h hs x sigma2 mu2
  
end % end of ROI
