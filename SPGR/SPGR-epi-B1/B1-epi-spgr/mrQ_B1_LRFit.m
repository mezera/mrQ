function mrQ_B1_LRFit(IDnum,jumpindex,jobindex)
%
% function mrQ_B1_LRFit(optName,jumpindex,jobindex)
%
%  This function, either through the SunGrid or alone, loads the relevant
%  data and fits the B1 in a specific region (volume). Herein, the imaging
%  volume region is also called the "box". The box is a location ranging in
%  size from a few hundred voxels to a few thousand.
%
%    ~INPUTS~
%         optName: This is optmization structure that was passed from
%                    mrQ_PD_LRB1SPGR_GridParams.m
%       jumpindex: How many boxes we fit in this grid call  (bookkeeping)
%        jobindex: The number of boxes it will start when calling the grid
%                  (bookkeeping)
%
%   ~OUTPUTS~
%  This function will save an output file with fitted parameters to a tmp
%  directory. This will be used later to make the B1 map.
%
% AM (C) Stanford University, VISTA
%

%% I. Initialization

[mrQpath]= mrQ_getPath(IDnum);
load(mrQpath);
optName=mrQ.B1.logname;

load(optName);
%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%check that this box has brain data
if ed>opt.N_Vox2Fit, ed=opt.N_Vox2Fit;end;

nIteration=ed-st+1;
%initialize parameters and saved outputs

%

% initialize the saved parameters
exitflag=zeros(nIteration,1);
resnorm=zeros(nIteration,1);
resnorm_noTresh=zeros(nIteration,1);

UseVoxN=zeros(nIteration,1);
skip=zeros(nIteration,1);
Iter=0;
B1=zeros(nIteration,1);
B1_noTresh=zeros(nIteration,1);
options = optimset('Display','off',...
    'MaxFunEvals',Inf,...
    'MaxIter',Inf,...
    'TolFun', 1e-6,...
    'TolX', 1e-10,...
    'Algorithm','trust-region-reflective'...
);
%    'Algorithm','levenberg-marquardt'

% make mask
BM=readFileNifti(opt.tisuuemaskFile);
BM=logical(BM.data);

% Get the needed parameters to fit
%EPI space
[tr, flipAngles,Res]=epiParams(opt);
SigMask=logical(ones(size(BM)));

N_Measure=length(tr);
ratios=nchoosek(1:N_Measure,2);


loc=find(SigMask);



%%  II. Go over box by box
tic    ;
for ii= st:ed,
    %run over the box you'd like to fit
    
    % clear parameters
    Iter= Iter+1;
    
    if ~(ii>length(loc))
        [S, t1, BM1,SZ, UseVoxN(Iter), skip(Iter), f ]=  mrQ_GetB1_LR_Data(opt,Res,BM,loc(ii));
    else
        skip(Iter)=1;
    end
    
    
    if  skip(Iter)==1
        %         disp(['skipping box bad data'])
    else
        % loop over Poly degrees
        S=reshape(S,prod(SZ(1:3)),SZ(4));
        
        f1=repmat(f(:),1,size(ratios,1));
        
        [B1(Iter) ,resnorm(Iter),~,exitflag(Iter)] = ...
            lsqnonlin(@(par) errB1_LR(par, flipAngles,tr,double(S),double(t1(:)),...
            double(f1),  BM1(:), ratios),...
            1,0.3,1.7,options);
        B1_noTresh(Iter)=B1(Iter);
        resnorm_noTresh(Iter)=resnorm(Iter);
        % Something wrong with the flip angle, I believe.
        % Check the ANat call and file orders!!!
        if  B1(Iter)==0.3 || B1(Iter)==0.7
            B1(Iter)=0;
            resnorm(Iter)=0;
        end
        
    end
    
    
end
toc
 %% III. Cross-Validation Fit
 % We can save some of the ratio and cross-validate the poly degree or box size.
       
name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'B1','B1_noTresh','resnorm','resnorm_noTresh','exitflag','st','ed','skip','UseVoxN' )

end

%%

function [tr, flipAngles,Res]=epiParams(opt)


% load infoormation
load (opt.AlignFile);

  % load information
    load (opt.AlignFile);

flipAngles=opt.FlipAngle;
tr=opt.TR;
end

%%