function logname=mrQ_PD_LRB1SPGR_GridParams(mrQ,FilterSize,pracent_coverage)
%%  logname=mrQ_PD_LRB1SPGR_GridParams(mrQ,AnalysisInfo,FilterSize,pracent_coverage)

% we building a stracture opt that will be used to fit B1   



opt.mrQfile=mrQ.name;

%opt.AnalysisInfoFile=AnalysisInfo.name;
%% get names and files 

outDir = mrQ.outDir;
subName=mrQ.sub;

%% LR sizeand mimumus local covarage
if notDefined('FilterSize')
    FilterSize =6;
   
end
if notDefined('pracent_coverage')
    pracent_coverage =0.33;
end

opt.pracent_coverage=pracent_coverage;
opt.FilterSize=FilterSize;



 opt.tisuuemaskFile=mrQ.maskepi_File;
 BM=readFileNifti(mrQ.maskepi_File);opt.pixdim=BM.pixdim;
 BM=logical(ones(size(BM.data)));
 opt.AlignFile=mrQ.spgr2epiAlignFile;



% use the mask where we will like to have a B1 mask
 opt.N_Vox2Fit=length(find(BM));




%%
sgename    = [subName '_B1'];
dirname    = [outDir '/tmpSGB1' ];
dirDatname = [outDir '/tmpSGB1dat'];
jumpindex  = 100; %number of boxs fro each SGR run

opt.dirDatname = dirDatname;
opt.name = [dirname '/B1boxfit_iter'] ;
opt.date = date;
opt.jumpindex = jumpindex;
opt.dirname=dirname;

opt.SGE=sgename;
% Save out a logfile with all the options used during processing
logname = [outDir '/fitLogB1.mat'];
opt.logname=logname;
%saving an information file we can load after if needed
save(opt.logname,'opt');


if clobber && (exist(dirname,'dir'))
    % in the case we start over and there are  old fits, so we will
    % deleat them
    eval(['! rm -r ' dirname]);
end