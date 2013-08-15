 outDir='/biac4/wandell/biac2/wandell2/data/WMDevo/phantom/spgr_32/20110812_0872/X_validationRegAnal'
 SunGrid=1;
 subName='phantom';
 proclass=1;
[opt]=mrQ_PD_multicoil_RgXv_GridCall(outDir,SunGrid,proclass);


    save(opt.logname,'opt');
clear
load /biac4/wandell/biac2/wandell2/data/WMDevo/phantom/spgr_32/20110812_0872/X_validationRegAnal/fitLog.mat
mrQ_CoilPD_gridFit(opt,opt.jumpindex,10); 

mrQ_CoilPD_gridFit(opt,opt.jumpindex,300);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,1848);

mrQ_CoilPD_gridFit(opt,opt.jumpindex,725);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,585);


%%
%%
clear
Strfilename='/biac4/wandell/biac2/wandell2/data/WMDevo/phantom/spgr_32/20110812_0872/X_validationRegAnal/fitLog.mat'
mrQ_buildPD(Strfilename)


%load('/biac4/wandell/biac2/wandell2/data/WMDevo/phantom/spgr_32/20110812_0872/X_validationRegAnal/TmpBoxRes.mat')
% %get the sizes
% A=sum(LinScaleMat,2);
% B=sum(LinScaleMat,1);
% wh1=find(A); wh2=find(B);
% Mat=LinScaleMat(wh1,wh2);
% % add a constant in the end
% Mat(end+1,1)=Mat(1,1)/(length(find(Mat(1,:)))-1);


%  y=zeros(size(Mat,1),1);
% y(end)=Mat(end,1);
% % %solve it as multi  linear eqation
% 
%  C=pinv(Mat'*Mat)*Mat'*y;
% %[U, ~] = eig(polyRatioMat'*polyRatioMat);
% C=pinv(Mat'*Mat)*Mat'*y;
% 
%  Mats=Mat(1:10,:);
% Mats(11,1)=1;
%  y=zeros(size(Mats,1),1);
% y(end)=Mats(end,1);
%  C=pinv(Mats'*Mats)*Mats'*y;

 %%  Brain
 clear
outDir='/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1'
 SunGrid=1;
 subName='AM';
 proclass=1;
[opt]=mrQ_PD_multicoil_RgXv_GridCall(outDir,SunGrid,proclass);


load /biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1/fitLog.mat
%mrQ_fitM0boxesCall(opt.logname,SunGrid,proclass);

mrQ_CoilPD_gridFit(opt,opt.jumpindex,300);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,183);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,205);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,1785);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,1488);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,1489);

%
clear
Strfilename='/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1/fitLog.mat'
mrQ_buildPD(Strfilename)



%%