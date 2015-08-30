function opt_Log_name=mrQ_buildPD(opt_Log_name,csffile,segfile,...
                                      RepErrThreshold,PrcCutOff,ErrorThresh)
% 
% opt=mrQ_buildPD_ver2(opt_Log_name,csffile,segfile,...
%                               RepErrThreshold,PrcCutOff,ErrorThresh)
%
% This function collects all the fitted coil gains in different small
% volumes (boxes) which were fitted in parallel, joins them back to a PD
% map, and calculates the WF map. 
%
% The fit can be done via one of three methods:
%  PDfit_Method =
%               = 1 : Use only T1 single-coil data to fit (default)
%               = 2 : Use only multi-coils to fit
%               = 3 : Use both T1 single-coil data *and* multi-coils to fit
%
% Building the WF map is done is 6 steps: 
%     Step 0: Loop and load all the fitted data of the different boxes, and 
%             check which of them have data.
%     Step 1: Get the x, y, z, and PD values of the boxes. 
%             --> mrQ_CalBoxPD_step1a.m
%     Step 2: Find the scalar between each box's PD.     
%             --> mrQ_ScaleBoxes_step2.m 
%     Step 3: Join the boxes to a PD image. 
%             --> mrQ_BoxJoinBox.m 
%     Step 4: Smoothe the gain values and re-calculate the PD. 
%             --> mrQ_smoothGain_step4b.m 
%     Step 5: Scale the PD CSF value to 1 to calculate the WF. This step 
%             uses the segmentation files csffile and segfile (inputs).
%             --> mrQ_PD2WF_step5.m
%
%
%   ~INPUTS~
%     opt_Log_name:   The name of the location of the "opt" structure
%          csffile:   The CSF segmentation file. If no csffile is selected,
%                              the default is to take from opt's outDir, 
%                              csf_seg_T1.nii.gz
%          segfile:   The T1-weighted tissue segmentation file. If no
%                              segfile is selected, the default is to take 
%                              from opt's outDir,T1w_tissue.nii.gz            
%  RepErrThreshold:   The acceptable amount of error when performing the
%                               multi-coils fit/refit of the data. It is a
%                               number between 0 and 1. Not applicable when
%                               PDfit_Method=1. Defaults to 0.6 when
%                               PDfit_Method=2 and to 0.1 when it is 3.
%        PrcCutOff:   The lower percentile cutoff for the regularization 
%                               weight [Default is 1]. Only for 
%                               PDfit_Method=3.
%      ErrorThresh:   The acceptable error when assessing the overlap when 
%                               pairing adjacent boxes. It is a number 
%                               between 0 and 1. Defaults to 0.05 when
%                               PDfit_Method=1 or 3, or to 0.01 when
%                               PDfit_Method=2 (values selected from
%                               previous experience with data)
%
%   ~OUTPUTS~
%              opt:   The updated opt structure of optimized parameters
%
% See also: mrQ_PD_multicoil_RgXv_GridCall.m 
%           mrQ_CoilPD_gridFit.m
%
% AM (C) Stanford University, VISTA
%
%

%% I. Load the fit information

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

if notDefined('ErrorThresh')
    if opt.PDfit_Method==1 || opt.PDfit_Method==3 % for T1 regularization
        ErrorThresh=0.01;
    else  %by simulation we expect more error and we must allow noisier fit
        ErrorThresh=0.05;
    end 
end

if notDefined('RepErrThreshold') %the repetition between coils max error
    if opt.PDfit_Method==3     % for  multi coils data with (3) or without (2) T1 regularization
        RepErrThreshold=0.1;
    elseif opt.PDfit_Method==2 % for Correlation regularization. By simulation, we expect more bias than in T1. We must allow it or we won't have data to work with.
        RepErrThreshold=0.6;%  multi coils
    else
        RepErrThreshold=1;
    end
end

if notDefined('PrcCutOff') %the lower percentile cutoff for the regularization weight
    PrcCutOff=1;
end

%% II. Step 0
% Get the fit for each box

% Loop over the fitted box files and load the parameters to matrices

jumpindex=opt.jumpindex; % number of box in each file
jobindexs=1:ceil(length(opt.wh)/jumpindex); % the numbers of saved files, each containing jumpindex boxes

if opt.PDfit_Method~=2 && opt.PDfit_Method~=1 % this output (cross-validation) is not part of the correlation method
   
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

% Select the boxes to combine

BoxesToUse=find(GoodBoxs)';

LowRepitionErr=BoxResidErr<RepErrThreshold;
BoxesToUse=find(GoodBoxs & LowRepitionErr)';

%figure;hist(BoxResidErr(BoxesToUse),1000)
% if opt.PDfit_Method~=2 && opt.PDfit_Method~=5
if opt.PDfit_Method~=2 && opt.PDfit_Method~=1

    RegW=BoxBestReg(BoxesToUse,2); LawW=prctile(RegW,PrcCutOff);
    BoxesToUse=find(GoodBoxs & LowRepitionErr &  BoxBestReg(:,2)>LawW)';
    % figure;plot(BoxBestReg(BoxesToUse),BoxResidErr(BoxesToUse),'*')
end

tmpfile=fullfile(opt.outDir,'Boxtmp');

%% III. Step 1
% Get the location and PD information for each box we will use

    [Boxes, PositiveBoxs,UnCorBoxs, UnSTDBoxs]=mrQ_CalBoxPD_step1a(opt,BoxesToUse,CoilGains);
    BoxesToUse=find(PositiveBoxs & UnSTDBoxs==0)';

%% IV. Step 2
% Find the ratio for each box's PD with its neighbors

[Boxes, ScaleMat]=mrQ_ScaleBoxes_step2(Boxes,BoxesToUse,opt,ErrorThresh);

% Solving the box weights by a system of linear equations that adjusts the
% median ratio between all the boxes
[Cbox SHub]=mrQ_boxScaleGlobLinear(ScaleMat);

%% V. Step 3
% Join the boxes to a PD image. 

[PD_fit,opt]= mrQ_BoxJoinBox(Boxes,Cbox,opt);
%save(tmpfile,'BoxesToUse','CoilGains','Boxes','ScaleMat','PD_fit')

%% VI. Step 4
% Get a smooth coil sensitivity in all locations, bring back to original image space, and calculate PD
[opt]=mrQ_smoothGain_step4b(opt,PD_fit); 

%% VII. Step 5
%Calculate WF by normalizing PD to CSF=1


%% calculate WF by normalize PD to CSF==1
% we are performing a different calibration now, carried out in mrQ_WF

% if unique( ~csffile==0)
%     opt=mrQ_PD2WF_step5(opt,csffile,segfile);
%     % I think it will only go in here if there are NO zeros in 'csffile' because [...?]
% end


