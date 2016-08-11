function ok = mrQ_PolyFitOrder(mask,data,flegUser)

if ~exist('flegUser','var')||isempty(flegUser),
    flegUser=0;
end;

addpath(genpath('/home/avivm/matlab/vistasoft/trunk/kendrick/kendrick'))
%addpath('~knk/code');

%load /home/avivm/matlab/hydrationLayer/trunk/ahpddata.mat

degs = 0:10;
params = {}; gains = {}; rs = {};
for p=1:length(degs)
 [params{p},gains{p},rs{p}] = fit3dpolynomialmodel(data,mask==1,degs(p));
end


if flegUser==100;
figure(10);
plot(degs,cat(2,rs{:}),'ro-');
xlabel('degree');
ylabel('correlation with the data');

p = input('select the polynomial degree ') ;
p=p+1; %degree 0:10 deg(p) are 1:11
data0 = zeros(size(data));
%data0(mask==1) = data(mask==1);

%for p=2:length(degs)

 basis = constructpolynomialmatrix3d(size(data),find(ones(size(data))),degs(p));
 ok = reshape(basis*params{p}',size(data));
 %drawnow; figure; imagesc(makeimagestack(data0)); axis equal tight;
 drawnow; figure(11); imagesc(makeimagestack(ok)); axis equal tight;
elseif flegUser~=100
        p=flegUser;
    else
  p=4;
  display=['using polynomial degree 4']
end;

  basis = constructpolynomialmatrix3d(size(data),find(ones(size(data))),degs(p));
 ok = reshape(basis*params{p}',size(data));

