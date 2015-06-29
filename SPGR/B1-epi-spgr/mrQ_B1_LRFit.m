function mrQ_B1_LRFit(opt,jumpindex,jobindex)
%
% mrQ_B1_gridFit(opt,jumpindex,jobindex)
%  this function call by the sun grid it load the relavant data and fit the
%  B1 in a rigion (volume).
% the imaging voulume region also call here "box". The box is a location (few voxel 100's to
% 1000's).  
% INPUTS:
%       opt - this is optmization structure that was passed from
%       jumpindex - how many boxes this grid call we fit (book keeping)
%       jobindex  - the number of box it will start whe ncalling the grid
%       (book keeping)
%
% OUTPUTS:
%  save an output file with fitted parameters in a tmp directorry
%   this will be used later  to make the B1 map

% SEE ALSO:
% AM (C) Stanford University, VISTA
%
%

%% I. Initialization



%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%cheack that this box have brain data
if ed>opt.N_Vox2Fit, ed=opt.N_Vox2Fit;end;

nIteration=ed-st+1;
%intilazie parameters and saved outputs

%


% initite the saved parameters
exitflag=zeros(nIteration,1);
resnorm=zeros(nIteration,1);
UseVoxN=zeros(nIteration,1);
skip=zeros(nIteration,1);
Iter=0;
B1=zeros(nIteration,1);

    options = optimset('Display','off',...
        'MaxFunEvals',Inf,...
        'MaxIter',Inf,...
        'TolFun', 1e-6,...
        'TolX', 1e-10,...
        'Algorithm','levenberg-marquardt');
  % make mask
BM=readFileNifti(opt.tisuuemaskFile); 
BM=logical(BM.data);    
  
    %% get the needed parameters to fit
    %EPI space
    [tr, flipAngles,Res]=epiParams(opt);
    SigMask=logical(ones(size(BM)));
    
    N_Measure=length(tr);
    ratios=nchoosek(1:N_Measure,2);



loc=find(SigMask);
%%  II. go over the box the boxs
       
for ii= st:ed,
    %run over the box you like to fit
   % clear parameters
    Iter= Iter+1;
    tic
   
    if ~(ii>length(loc))
          [S, t1, BM1,SZ, UseVoxN(Iter), skip(Iter), f ]=  mrQ_GetB1_LR_Data(opt,Res,BM,loc(ii));
    else
        skip(Iter)=1;
    end
      
          
    
    
    
    if  skip(Iter)==1
        disp(['skipping box  bad data'])
         
    else
        %% lop over Poly degrees
     S=reshape(S,prod(SZ(1:3)),SZ(4));
            
        f1=repmat(f(:),1,size(ratios,1));
        
  [B1(Iter) ,resnorm(Iter),~,exitflag(Iter)] = ...
    lsqnonlin(@(par) errB1_LR(par, flipAngles,tr,double(S),double(t1(:)),...
     double(f1),  BM1(:), ratios),...
    1,[],[],options);
        
%% some thing wrong with the flip angle i bellive. cheack the ANat call. and file orders!!!
 
     end
     
toc
        end; 
        %%  X-Validation Fit we can save some of the ratio and corrsvalidate the poly degree or  box size .
       
        %




name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'B1','resnorm','exitflag','st','ed','skip','UseVoxN' )

end

%%

function [tr, flipAngles,Res]=epiParams(opt)

  % load infoormation
    load (opt.AlignFile);
flipAngles=ResInfo.flipAngles;
tr=ResInfo.tr;
end


%%



