function [f,k,s0,gof] = mrQ_relaxFitMtDist(data,delta,t1,s0,tr,flipAngle,brainMask,fitMethod,outMontage)
%
% [f,k,s0,gof] = relaxFitMtDist(data,delta,t1,s0,tr,flipAngle,brainMask,fitMethod,outMontage)
% 
% Computes a nonlinear fit of the bound-pool (f) map.
%
% Returns:
%
% SEE ALSO:
% 
% relaxFitT1.m to fit the t1 and pd maps used by this function.
%
% HISTORY:
% 2008.02.26 RFD: wrote it.

if(~exist('fitMethod','var')||isempty(fitMethod))
    fitMethod = 'fmin';
end


opt.flipAngle = flipAngle*pi/180;
opt.tr = tr/1000;
opt.sz = size(data);


brainInds = find(brainMask);

%opt.numVoxelsPerUpdate = 16000;

%opt.nVoxAll = length(brainInds);

 for(ii=1:size(data,4))
   tmpVol = data(:,:,:,ii);
   tmpMT(ii,:) = tmpVol(brainInds);
 end
 clear tmpVol;
 data = tmpMT;
 clear tmpMT;

opt.r1 = 1./t1(brainInds);
opt.s0 = s0(brainInds);
opt.wh   = find(brainMask);
opt.B1   = B1(brainMask);

opt.outDir = [outDir '/tmpSG'];
opt.lb     = [0 0];
opt.ub     = [inf 10000];
opt.name   = '/BPF_fit';
% f = zeros(1,nVoxAll); 
% k = zeros(1,nVoxAll);
% gof = zeros(1,nVoxAll);
% totalSecs = 0;

%tmpImg = zeros(sz(1:3));
%tmpName = tempname
%nSteps = ceil(nVoxAll/numVoxelsPerUpdate);
%totalTime = 0;
%tic;

if SGE==1;
    jumpindex=1000;
    if (~exist([outDir '/tmpSG'],'dir')), mkdir([outDir '/tmpSG']);
        % the result form the grid will be saved in a tmporery directory
        sgerun('mrQ_relaxFitMt.m(opt,1000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        
    else









for(ii=1:nSteps)
    if(useParfor), matlabpool; end
    curInd = (ii-1)*numVoxelsPerUpdate+1;
    endInd = min(curInd+numVoxelsPerUpdate,nVoxAll);
    [tf,tk,tgof] = relaxFitMt(data(:,curInd:endInd),delta,r1(curInd:endInd),s0(curInd:endInd),tr,flipAngle,fitMethod,useParfor);
    prevSecs = toc;
    totalTime = totalTime+prevSecs;
    f(curInd:endInd) = tf;
    k(curInd:endInd) = tk;
    gof(curInd:endInd) = tgof;
    secsPerVox = prevSecs/(endInd-curInd+1);
    estTime = secsPerVox*(nVoxAll-endInd);
    if(estTime>5400) estTime=estTime./3600; estTimeUnits='hours';
    elseif(estTime>90) estTime=estTime./60; estTimeUnits='minutes';
    else estTimeUnits='seconds'; 
    end
    fprintf('Processed %d of %d voxels- %0.1f %s remaining (%0.3f secs per vox)...\n',endInd,nVoxAll,estTime,estTimeUnits,secsPerVox);
    tmpImg(brainInds) = f;
    m = makeMontage(tmpImg);
    if(max(f)>0), m = m./max(f); end
    m = uint8(round(m.*255));
    imwrite(m,outMontage);
    save(tmpName,'k','f','gof','brainMask');
    tic;
    if(useParfor), matlabpool close; end
end
fprintf('Processed %d voxels in %0.2f hours.\n',totalTime/3600);

% f = vertcat(tf);
% k = vertcat(tk);
% gof = vertcat(tgof);

im=zeros(size(brainMask)); im(brainInds) = f; f = im;
im=zeros(size(brainMask)); im(brainInds) = k; k = im;
im=zeros(size(brainMask)); im(brainInds) = gof; gof = im;
im=zeros(size(brainMask)); im(brainInds) = s0; s0 = im;

return;







opt.wh   = find(brainMask);
opt.Gain = Gain(brainMask);
opt.B1   = B1(brainMask);

opt.outDir = [outDir '/tmpSG'];
opt.lb     = [0 0];
opt.ub     = [inf 10000];
opt.name   = '/T1PDlsqVx';


%% Perform the optimization (optionally using the SGE)

% USE THE SGE
if SGE==1;
    jumpindex=2000;
    if (~exist([outDir '/tmpSG'],'dir')), mkdir([outDir '/tmpSG']);
        % the result form the grid will be saved in a tmporery directory
        sgerun('mrQ_fitT1PD_SGE(opt,2000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        
    else
        an1 = input( 'Unfinished SGE run found: Would you like to try and finish the existing sge run? Press 1 if yes. To start over press 0 ');
        
        % Continue existing SGE run from where we left it last time
        % we find the fit that are missing
        if an1==1
            reval=[];
            list=ls(opt.outDir);
            ch=[1:jumpindex:length(opt.wh)];
            k=0;
            for ii=1:length(ch),
                
                ex=['_' num2str(ch(ii)) '_'];
                if length(regexp(list, ex))==0,
                    k=k+1;
                    reval(k)=(ii);
                end
            end
            
            if length(find(reval))>0
                % clean the sge output dir and run the missing fit
                eval(['!rm -f ~/sgeoutput/*' sgename '*'])
                sgerun('mrQ_fitT1PD_SGE(opt,2000,jobindex);',sgename,1,reval);
                
            end
            list=ls(opt.outDir);
            
            % Restart the SGE processing from the beginning
        elseif an1==0
            t=pwd;
            cd (outDir)
            !rm -rf tmpSG
            cd (t);
            
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
            mkdir([outDir '/tmpSG']);
            sgerun('mrQ_fitT1PD_SGE(opt,2000,jobindex);',sgename,1,1:ceil(length(opt.wh)/jumpindex));
        else
            error;
        end
    end
    
    %% build the data that was fit by the SGE to a T1 nd M0 maps
    % This loop checks if all the outputs have been saved and waits until
    % they are all done
    StopAndSave=0;
    fNum=ceil(length(opt.wh)/jumpindex);
    tic
    while StopAndSave==0
        % List all the files that have been created from the call to the
        % grid
        
        list=ls(opt.outDir);
        % Check if all the files have been made.  If they are, then collect
        % all the nodes and move on.
        if length(regexp(list, '.mat'))==fNum,
            StopAndSave=1;
            
            % Loop over the nodes and collect the output
            for i=1:fNum
                st=1 +(i-1)*jumpindex;
                ed=st+jumpindex-1;
                
                if ed>length(opt.wh), ed=length(opt.wh);end
                
                name=[opt.outDir '/' opt.name '_' num2str(st) '_' num2str(ed) '.mat'];
                load (name);
                t11(st:ed)=res(2,:);
                pd1(st:ed)=res(1,:);
                resnorm1(st:ed)=resnorm;
                
            end
            % Once we have collected all the nodes we delete the temporary
            t=pwd;
            cd (outDir)
            !rm -r tmpSG
            cd (t);
            eval(['!rm -f ~/sgeoutput/*' sgename '*'])
        end
        
        
        % Record how much time has elapsed since the call to the grid.
        t = toc;
        % If too much time has elapsed then we recall the grid;
        if t > 86400% 24hours
            reval=[]
            ch=[1:jumpindex:length(opt.wh)]; %the nude filre name
            k=0;
            reval=[];
            
            for ii=1:length(ch),
                
                ex=['_' num2str(ch(ii)) '_'];
                if length(regexp(list, ex))==0,
                    k=k+1;
                    reval(k)=(ii); % we make a list of the grid run that are not done yet
                    
                end
            end;
            if length(find(reval))>0
                eval(['!rm ~/sgeoutput/*' sgename '*']) % we delete all our relevant grid jobs
                sgerun('mrQ_fitT1PD_SGE(opt,500,jobindex);',sgename,1,reval,[],[],3000); % we run the missing oupput again
            else
                display('somting is wrong in SGE run')
                error
            end
        end
        
    end
    
    % NO SGE
    %using the local computer to fit T1 and the sunGrid
else
    % Run the optimization without using the SGE
    for i= 1:length(opt.wh)
        
        [res(:,i), resnorm(i)] = lsqnonlin(@(par) errT1PD(par,opt.flipAngles,opt.tr,opt.s(i,:),opt.Gain(i),opt.B1(i),1,[]),opt.x0(i,:),opt.lb,opt.ub,options);
        
    end
    
    t11(:)=res(:,2);
    pd1(st:ed)=res(:,1);
    
end

T1 = zeros(size(brainMask));
PD = T1; resNorm=PD;
T1(opt.wh) = t11(:)./1000;
PD(opt.wh) = pd1(:);
resNorm(opt.wh) = resnorm1(:);

