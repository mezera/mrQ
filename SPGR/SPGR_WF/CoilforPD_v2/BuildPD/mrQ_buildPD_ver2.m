function opt_Log_name=mrQ_buildPD_ver2(opt_Log_name,csffile,segfile,RepErrTreshold,PrcCutOff,ErrorTresh)
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
    segfile = fullfile(opt.outDir, 'T1w_tissue.nii.gz');
end

if notDefined('ErrorTresh')
    if opt.PDfit_Method==1 % for T1 regularization
        ErrorTresh=0.01;
    else  %by simulation we expect more error and we must allow noisyer fit
        ErrorTresh=0.05;
    end
    
end

if notDefined('RepErrTreshold') %the repetition betwwen coils max error
    if opt.PDfit_Method==1 % for T1 regularization
        
        RepErrTreshold=0.1;
    elseif opt.PDfit_Method==2 % for Correlation regularization by simulation we expect more bias then T1. we must allow it or we won't have data to work with
        RepErrTreshold=0.6;
    elseif opt.PDfit_Method==3 % for ridge regularization by simulation we expect very strong bias.we must allow it or we won't have data to work with
        RepErrTreshold=0.95;
    else
        RepErrTreshold=1;
    end
end

if notDefined('PrcCutOff') %the lower pracential cut off for the regularization wight
    PrcCutOff=1;
end
%% get the fit for each box
% loop over the fitted box files and load the parameters to  matrixs
jumpindex=opt.jumpindex; % number of box in each file
jobindexs=1:ceil(length(opt.wh)/jumpindex);  % the numbers of saved file each with jumpindex boxs

% if opt.PDfit_Method~=2 & opt.PDfit_Method~=5 % this oupput (X-vlidation ) are not part of the correlation method
if opt.PDfit_Method~=2 & opt.PDfit_Method~=1 % this output (X-vlidation ) are not part of the correlation method
   
GoodBoxs=zeros(length(opt.wh),1);
       BoxBestReg=zeros(length(opt.wh),2);
    BoxExitflag=zeros(length(opt.wh),1);
    BoxX_valdationErr=zeros(2,length(opt.lambda),length(opt.wh));
    BoxX_ValSign=zeros(length(opt.wh),1);
    BoxResidErr=zeros(length(opt.wh),1);
    
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
                BoxBestReg(boxes(jj),:)=opt.lambda(BestReg(jj,:));
                BoxX_ValSign(boxes(jj))=X_valdationErrSN(jj);
                BoxX_valdationErr(:,:,boxes(jj))=X_valdationErr(:,:,jj);
                BoxExitflag(boxes(jj))=exitflag(jj);
                coils=Clists(:,jj);
                coils=coils(coils~=0);
                Nc=length(coils);
                CoilGains(boxes(jj)).Clist=coils;
                CoilGains(boxes(jj)).g=gEst(:,1:Nc,jj);
                BoxResidErr(boxes(jj))=ResidErr(jj,:);
            end
        end
    end
    
end

if opt.PDfit_Method==2
    GoodBoxs=zeros(length(opt.wh),1);
    BoxExitflag=zeros(length(opt.wh),1);
    BoxResidErr=zeros(length(opt.wh),1);
    
    for ii=1:length(jobindexs);
        
        % the number of the box in the file (also part of the saved file name)
        st=1 +(jobindexs(ii)-1)*jumpindex;
        ed=st+jumpindex-1;
        if ed>length(opt.wh), ed=length(opt.wh);end;
        
        
        %% this is temporary fix !!!
        
        FileName=[ opt.name '_' num2str(st) '_' num2str(ed)];
        
        load(FileName);
        boxes=[st:ed];
        
        for jj=1:length(boxes)
            if skip(jj)==0
                GoodBoxs(boxes(jj))=1;
                BoxExitflag(boxes(jj))=exitflag(jj);
                coils=Clists(:,jj);
                coils=coils(coils~=0);
                Nc=length(coils);
                CoilGains(boxes(jj)).Clist=coils;
                CoilGains(boxes(jj)).g=gEst(:,1:Nc,jj);
                BoxResidErr(boxes(jj))=ResidErr(jj,:);
            end
        end
    end
    
end
%%

% if opt.PDfit_Method==5
if opt.PDfit_Method==1
    GoodBoxs=zeros(length(opt.wh),1);
    BoxResidErr=zeros(length(opt.wh),1);
    
    
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
                
                CoilGains(boxes(jj)).Clist=1;
                CoilGains(boxes(jj)).g=gEst(:,jj);
                BoxResidErr(boxes(jj))=-100;
            end
        end
    end
end


%%

%% select the box to combine

BoxesToUse=find(GoodBoxs)';


LowRepitionErr=BoxResidErr<RepErrTreshold;
BoxesToUse=find(GoodBoxs & LowRepitionErr)';
%figure;hist(BoxResidErr(BoxesToUse),1000)
% if opt.PDfit_Method~=2 && opt.PDfit_Method~=5
if opt.PDfit_Method~=2 && opt.PDfit_Method~=1

    RegW=BoxBestReg(BoxesToUse,2); LawW=prctile(RegW,PrcCutOff);
    BoxesToUse=find(GoodBoxs & LowRepitionErr &  BoxBestReg(:,2)>LawW)';
    % figure;plot(BoxBestReg(BoxesToUse),BoxResidErr(BoxesToUse),'*')
end



tmpfile=fullfile(opt.outDir,'Boxtmp');
%save(tmpfile,'BoxesToUse','CoilGains')
%% get the location and PD information to each box we will use
if opt.PDfit_Method==1
    [Boxes, PositiveBoxs,UnCorBoxs, UnSTDBoxs]=mrQ_CalBoxPD_step1a(opt,BoxesToUse,CoilGains);
    BoxesToUse=find(PositiveBoxs & UnSTDBoxs==0)';
    %[Boxes, PositiveBoxs]=mrQ_CalBoxPD_step1(opt,BoxesToUse,CoilGains);
    %Boxe5sToUse=find(PositiveBoxs)';
elseif opt.PDfit_Method~=1
    [Boxes, PositiveBoxs,UnCorBoxs, UnSTDBoxs]=mrQ_CalBoxPD_step1a(opt,BoxesToUse,CoilGains);
    BoxesToUse=find(PositiveBoxs & UnSTDBoxs==0)';
    %this not use the Corr error anyway, right???
    %BoxesToUse=find(PositiveBoxs & UnCorBoxs==0 & UnSTDBoxs==0)';
    %save(tmpfile,'BoxesToUse','CoilGains','Boxes')
end
%%  join the boxs
% find the ratio for each boxes PD with it's niebores

[Boxes, ScaleMat]=mrQ_ScaleBoxes_step2(Boxes,BoxesToUse,opt,ErrorTresh);
%save('tmp','LinScaleMat','ScaleMat')
%save(tmpfile,'BoxesToUse','CoilGains','Boxes','ScaleMat')

% solving the box wights by system of linear eqation that adjust the median ratio
% between all the boxs
[Cbox SHub]=mrQ_boxScaleGlobLinear(ScaleMat);

% join the box acording to the box Costat calculate by
% mrQ_boxScaleGlobLinear.m
[PD_fit,opt]= mrQ_BoxJoinBox(Boxes,Cbox,opt);

%Try to use this function if is the other didn't work. This will join the boxes one by one and reject outlayer.
%this function may be dependent on the starting box.
%[PD_fit]= mrQ_BoxsWiseAlignment_step3(Boxes,ScaleMat,BoxesToUse,BoxesToUse(1),opt);
save(tmpfile,'BoxesToUse','CoilGains','Boxes','ScaleMat','PD_fit')

%% get a smooth G and in all location, bring back to original image space anscalculate P .
[opt]=mrQ_smoothGain_step4b(opt,PD_fit); 
%eval(['! rm ' tmpfile])

%% calculate WF by normalize PD to CSF==1
if unique( ~csffile==0) %% i think it will only go in here if ther are NO zeros in 'csffile'because 
    opt=mrQ_PD2WF_step5_ver2(opt,csffile,segfile);
end


