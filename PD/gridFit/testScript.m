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


%%
outDir='/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1'
 SunGrid=1;
 subName='AM';
 proclass=1;
[opt]=mrQ_PD_multicoil_RgXv_GridCall(outDir,SunGrid,proclass);

clear
load /biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1/fitLog.mat

mrQ_CoilPD_gridFit(opt,opt.jumpindex,300);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,183);
mrQ_CoilPD_gridFit(opt,opt.jumpindex,205);


%%