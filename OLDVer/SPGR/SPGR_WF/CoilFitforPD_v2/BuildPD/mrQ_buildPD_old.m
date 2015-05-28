function opt=mrQ_buildPD_old(opt_Log_name,csffile)
% this function colect all the fitted coil gains in different small voulumes (box) that was fitted in paraller
%and join them back to a PD map and calculate the WF map. Strfilenameis the log file of the parralel fits.
%   opt=mrQ_buildPD(opt_Log_name,csffile,segfile)
% building the WF is done is 6 steps:
% step 0: loop and load all the fitted data of the different boses. and cheack which of them have data
% step 1- get the x,y,z and pd values of the boxes.        --> mrQ_CalBoxPD_step1.m
% step 2 - % find the scaler btween each boxes PD      --> mrQ_ScaleBoxes_step2.m
% step 3 - join the boxes to an PD image.                     --> mrQ_BoxsWiseAlignment_step3.m
% step 4 smooth the Gain values and re-calculate PD --> mrQ_smoothGain_step4.m
% step 5 scale PD csf value to 1 to calculate the WF.   --> mrQ_PD2WF_step5.m
%            step 5 use segmentation files: csffile,segfile. 
%
% See also mrQ_PD_multicoil_RgXv_GridCall.m & mrQ_CoilPD_gridFit.m 
%
% AM (C) Stanford University, VISTA
%



%% load  the fit information
if (~exist(opt_Log_name,'file')),
    
    disp(['cant find file : ' opt_Log_name  ])
    error
else
    load  (opt_Log_name)
end

if notDefined('csffile')
    csffile = fullfile(opt.outDir, 'csf_seg_T1.nii.gz');
end
if notDefined('segfile')
    segfile = fullfile(opt.outDir, 't1_bet_seg.nii.gz');
end

%% get the fit for each box
% loop over the fitted box files and load the parameters to  matrixs
jumpindex=opt.jumpindex; % number of box in each file
jobindexs=1:ceil(length(opt.wh)/jumpindex);  % the numbers of saved file each with jumpindex boxs

GoodBoxs=zeros(length(opt.wh),1);
BoxBestReg=zeros(length(opt.wh),2);
BoxExitflag=zeros(length(opt.wh),1);
BoxX_valdationErr=zeros(2,length(opt.lambda),length(opt.wh));
BoxX_ValSign=zeros(length(opt.wh),1);


for ii=1:length(jobindexs);
    
    % the number of the box in the file (also part of the saved file name)
    st=1 +(jobindexs(ii)-1)*jumpindex;
    ed=st+jumpindex-1;
    if ed>length(opt.wh), ed=length(opt.wh);end;
    
    FileName=[ opt.name '_' num2str(st) '_' num2str(ed)];
    
    load(FileName);
    boxes=[st:ed];
    
    for jj=1:length(boxes)
        if skip(jj)==0
            GoodBoxs(boxes(jj))=1;
            BoxBestReg(boxes(jj),:)=BestReg(jj,:);
            BoxX_ValSign(boxes(jj))=X_valdationErrSN(jj);
            BoxX_valdationErr(:,:,boxes(jj))=X_valdationErr(:,:,jj);
            BoxExitflag(boxes(jj))=exitflag(jj);
            coils=Clists(:,jj);
            coils=coils(coils~=0);
            Nc=length(coils);
            CoilGains(boxes(jj)).Clist=coils;
            CoilGains(boxes(jj)).g=gEst(:,1:Nc,jj);
        end
    end
end


%%
BoxMaxReg=zeros(length(opt.wh),2);
for ii=1:length(opt.wh)
    
    
    if GoodBoxs(ii)==1;
        X_valdationErr=BoxX_valdationErr(:,:,ii);
        for jj=1:2
            [v, ind]=sort(   X_valdationErr(jj,:)./min(X_valdationErr(jj,:)));
            
            BoxMaxReg(ii,jj)=min(ind( find(v<1.05)));
        end
    end
end
BoxesToUse=find(GoodBoxs)';

tmpfile=fullfile(opt.outDir,'Boxtmp');
%save(tmpfile,'BoxesToUse','CoilGains')
%% get the location and PD information to each box we will use

[Boxes, PositiveBoxs]=mrQ_CalBoxPD_step1(opt,BoxesToUse,CoilGains);
BoxesToUse=find(PositiveBoxs)';

%save(tmpfile,'BoxesToUse','CoilGains','Boxes')

%%  join the boxs
% find the ratio for each boxes PD with it's niebores
[Boxes, ScaleMat]=mrQ_ScaleBoxes_step2(Boxes,BoxesToUse,opt);
%save('tmp','LinScaleMat','ScaleMat')
%save(tmpfile,'BoxesToUse','CoilGains','Boxes','ScaleMat')

% solving the box wights by system of linear eqation that adjust the median ratio
% between all the boxs 
[Cbox SHub]=mrQ_boxScaleGlobLinear(ScaleMat);

% join the box acording to the box Costat calculate by
% mrQ_boxScaleGlobLinear.m
[PD_fit]= mrQ_BoxJoinBox(Boxes,Cbox,opt);

%Try to use this function if is the other didn't work. This will join the boxes one by one and reject outlayer.
%this function may be dependent on the starting box.
%[PD_fit]= mrQ_BoxsWiseAlignment_step3(Boxes,ScaleMat,BoxesToUse,BoxesToUse(1),opt);
%save(tmpfile,'BoxesToUse','CoilGains','Boxes','ScaleMat','PD_fit')

%% get a smooth G and in all location, bring back to original image space anscalculate P . 
[opt]=mrQ_smoothGain_step4(opt,PD_fit);
%eval(['! rm ' tmpfile])

%% calculate WF by normalize PD to CSF==1
if unique( ~csffile==0)
opt=mrQ_PD2WF_step5(opt,csffile,segfile);
end