function T1FitDisplayBrain(loadStr, saveStr, method)

% T1FitDisplayBrain(loadStr, saveStr, method)
%  
%    loadStr: the data that the fit was performed on, to be loaded
%    saveStr: the results of the fit to be loaded
%    method:  what fitting method to use,
%             NLS or NLSPR
%
% A simple segmentation method is used to get WM and GM peak values.
% More elaborated methods could be used, but that was beyond the scope of
% our work. 
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University
% 

%%
return;


%% %%%%%% Nothing in here is run becuase of the return just above. %%%%%%%%
%% we don't not use this  set of plot and figure wany more

% Threshold for region growing; 300 at 1.5T, 400 at 3T
thresh = 400;

% Load the original data 
load(loadStr);

dims = size(data); 

% See how many slices there are
if numel(dims) > 3
  nbslice = dims(3); 
else
  nbslice = 1;
   % Make data a 4-D array regardless of number of slices
  tmpData(:,:,1,:) = data;
  data = tmpData;
  clear tmpData;
end

load(saveStr);         % This loads ll_T1
T1_original   = ll_T1; % Store it so ll_T1 can be a single slice;
Mask_original = mask;  % Store it so mask can be a single slice

% Inserting a short pause, otherwise some computers have problems
pause(1);

% Whole brain histogram 
T1roi = T1_original(:,:,:,1);
T1roi(find(T1roi==extra.T1Vec(1) | T1roi==extra.T1Vec(end))) = 0;

x  = min(T1roi(find(T1roi>0))):1:max(T1roi(find(T1roi>0)));
x  = x';
h  = histc(T1roi(find(T1roi>0)),x); % centered on integers
hs = medfilt1(h,10);                % 10 ms median filter cf. Tozer03 or Tofts

figure
bar(x,hs,1)
ylabel('Number of pixels')
xlabel('T_1 [ms]')
customFormat

[Nx,Ny] = size(T1roi); 


% Skull
roiname = 'Skull';
Lskull = ones(Nx,Ny,3); 

for zz = 1:nbslice
	ll_T1 = squeeze(T1_original(:,:,zz,:));
	mask = squeeze(Mask_original(:,:,zz));
	
	% Mask out points where grid search hit boundary
	T1 = ll_T1(:,:,1);
	mask(find(T1==extra.T1Vec(1) | T1==extra.T1Vec(end))) = 0;
	T1 = T1.*mask;
	
	satisfied = 0;
	while (satisfied == 0)
		disp(['Select seed in ' roiname ' for slice ' num2str(zz)])
		Lskull(:,:,zz) = regiongrowing(T1,thresh);
		figure
		imshow(Lskull(:,:,zz))
		satisfied = input([roiname ' properly selected? 1: yes, 0: no --- ']);
	end
	
end

% WM histogram to get peak value

roiname = 'Gray Matter';
T1roiWM = zeros(Nx,Ny,3);
T1roiGM = zeros(Nx,Ny,3);
L3 = ones(Nx,Ny,3);

for zz = 1:nbslice
	ll_T1 = squeeze(T1_original(:,:,zz,:));
	mask = squeeze(Mask_original(:,:,zz));
	
	% Mask out points where grid search hit boundary
	T1 = ll_T1(:,:,1);
	mask(find(T1==extra.T1Vec(1) | T1==extra.T1Vec(end))) = 0;
	T1 = T1.*mask;
	
	disp(['Select seed in ' roiname]) 
	J = regiongrowing(T1,thresh);
	L3(:,:,zz) = L3(:,:,zz)-J;
	figure
	imshow(J)

	title(roiname)
	T1roiWM(:,:,zz) = T1.*J;
end
T1roi = [T1roiWM(:,:,1),T1roiWM(:,:,2),T1roiWM(:,:,3)]; 
	
% Check if all elements are zero
if all( all(T1roi == 0) )
	error('All T1 estimates in the ROI is zero, no histogram can be plotted')
end
  
% Display a histogram of T1 in the ROI.  
% Bin width = Grid precision

x = min(T1roi(find(T1roi>0))):1:max(T1roi(find(T1roi>0)));
x = x';
h = histc(T1roi(find(T1roi>0)),x); % centered on integers
hs = medfilt1(h,10); % 10ms median filter cf. Tozer03 or Tofts
mu2 = mean(x(find(hs==max(hs))));
sigma2 = round(sqrt(1/sum(hs)*sum(hs.*((x - mu2).^2))));

figure
bar(x,hs,1)
titlename = [roiname, ': ', num2str(sum(h)), ' pixels, mode ', num2str(round(mu2)),...
    ' ms, \sigma = ', num2str(round(sigma2))]
title(titlename)
ylabel('Number of pixels')
xlabel('T_1 [ms]')
customFormat
clear h hs x sigma2 mu2

% GM histogram to get peak value
roiname = 'White Matter';

for zz = 1:nbslice
	ll_T1 = squeeze(T1_original(:,:,zz,:));
	mask = squeeze(Mask_original(:,:,zz));
	% Mask out points where grid search hit boundary
	T1 = ll_T1(:,:,1);
	mask(find(T1==extra.T1Vec(1) | T1==extra.T1Vec(end))) = 0;
	T1 = T1.*mask;

	L = L3(:,:,zz);
	L(find(Lskull(:,:,zz)==1)) = 0;
	L(find(T1==0))=0;
	figure, imshow(L,[])
	title(roiname)
	customFormat

	T1roiGM(:,:,zz) = T1.*L;
end
T1roi = [T1roiGM(:,:,1),T1roiGM(:,:,2),T1roiGM(:,:,3)]; 
  
% Check if all elements are zero
if all( all(T1roi == 0) )
	error('All T1 estimates in the ROI is zero, no histogram can be plotted')
end
  
% Display a histogram of T1 in the ROI.  
% Bin width = Grid precision
T1roi(find(T1roi>850|T1roi<550)) = 0; % empirical way to remove what is not WM 
x = min(T1roi(find(T1roi>0))):1:max(T1roi(find(T1roi>0)));
x = x';
h = histc(T1roi(find(T1roi>0)),x); % centered on integers
hs = medfilt1(h,10); % 10ms median filter cf. Tozer03 or Tofts
mu2 = mean(x(find(hs==max(hs))));
sigma2 = round(sqrt(1/sum(hs)*sum(hs.*((x - mu2).^2))));

figure
bar(x,hs,1)
titlename = [roiname, ': ', num2str(sum(h)), ' pixels, mode ', num2str(round(mu2)),...
    ' ms, \sigma = ', num2str(round(sigma2))]
title(titlename)
ylabel('Number of pixels')
xlabel('T_1 [ms]')
customFormat
clear h hs x sigma2 mu2

