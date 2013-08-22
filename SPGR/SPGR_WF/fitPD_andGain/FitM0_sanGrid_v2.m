function FitM0_sanGrid_v2(opt,jumpindex,jobindex)
%
% FitM0_sanGrid_v2(opt,jumpindex,jobindex)
%  this function call by the sun grid it load the relavant data and fit the
%  PD and coils bias of image rigion.
%the image region also call here box. is a location (few voxel 100's to
%1000's).  the idea is to use the information in man coils. we tery to find
%to PD that is similar to all the coils images and fit the coil bias that
%is diferent for each coil.
%  The fit is done in few step. 1. we get the box data.
%2. we fit subset of coils (nested function )
%3. we do qualtiy check to see if it worked o.k.  
%4. we save the resutl
%
% INPUTS:
%       opt - this is optmization structure that was passed from
%       mrQ_fitPD_multicoil and it have all the needed information
%       jumpindex - how many boxes this grid call we fit (book keeping)
%       jobindex  - the number of box it will start whe ncalling the grid
%       (book keeping)
%
% OUTPUTS:
%  save an output file with fitted parameters in a tmp directorry
%   this will be used lster by mrQfitPD_multiCoils_M0 to make the PD map

% SEE ALSO:
% mrQ_fitPD_multicoil
% (C) Stanford University, VISTA
%
%

%% I. Initialization

% use only one cpu or it will be very hard to work on the computer we used.
%maxNumCompThreads(1);
%maxnumcompthreads(1);

%find the box to work on
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;

%cheack that this box have brain data
if ed>length(opt{1}.wh), ed=length(opt{1}.wh);end;

%load the multi coil information
M=readFileNifti(opt{1}.dat);
coils=size(M.data,4);

%this is away to clear out layer but it can make problem so we can think
%about it again
Mall=mean(M.data,4);
...
    Mall(Mall<mean(Mall(opt{1}.brainMask))-3*std((Mall(opt{1}.brainMask))))=0;
Mall=logical(Mall);
for i=1:coils,
    tmp=M.data(:,:,:,i);
    tmp(~Mall)=0;
    M.data(:,:,:,i)=tmp;
end;

% make the polynomial that will be use to fit the coil gains
    [Poly,str] = constructpolynomialmatrix3d(opt{1}.boxS,find(ones(opt{1}.boxS)),opt{1}.degrees);
%[Poly,str] = constructpolynomialmatrix3d_1(opt{1}.boxS,find(ones(opt{1}.boxS)),[],[],opt{1}.str);

opt{1}.Poly=Poly;
opt{1}.lb = ones(opt{1}.numIn,size(opt{1}.Poly,2)).*-inf;
opt{1}.ub = ones(opt{1}.numIn,size(opt{1}.Poly,2)).*inf;

    %intilazie parameters and saved outputs
nn=1;
res=[];resnorm=0;
exitflag=-100;
coilList=0;skip=0;
coilList=zeros(opt{1}.numIn,ed-st+1);
    
%%  II. go over the box the boxs


for i= st:ed,   
     %run over the box you like to fit
  %%
    tic
    %find the x,y,z location of the box (this is not x,y,z location in image space but
    %grid of boxes we made by mashgrid in  mrQ_fitPD_multicoil
    [fb(i,1,1) fb(i,1,2) fb(i,1,3)]=ind2sub(size(opt{1}.X),opt{1}.wh(i));
    
    % get all the relevant box information raday to be fitted
    [x01,inDat1,bm1,szB,skip(nn) coils_list]= mrQM0fitini_Param_vv(opt{1},fb(i,:),M,coils);
    
    %cheack that there is usabale coils information if not we won't fit this box 
    if (size(inDat1,2))<4 
        skip(nn)=1;
        %no point to try; there is no data to even try to fit (less then 4 coil usefull data).
    end
    
    if skip(nn)==0
        % we will work will subset of the coils as define in opt (typicaly 8) in the casethere are more avialbe coils 
        use=[1:min(size(inDat1,2),opt{1}.numIn)];
        
        %Fit the data this with  a nested function (below)
        [res(:,:,nn) resnorm(nn) exitflag(nn) num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt);
        
        %check that the number of corect fitted coil are ther number that are in
        if num==opt{1}.numIn
            % if yes; all done we move on; save the coil we end up fitting
            coilList(1:length(use),nn)=coils_list(use);
        else%if num<5,
            %some coils are hard to fit so we it better replaced if we have
            %more coils then we used. or we can just move them and fit with
            %less
            
            replace=opt{1}.numIn-num;
            extra=size(inDat1,2)-opt{1}.numIn;
            next=opt{1}.numIn+1;
            %check that there are exstra coils 
            while (replace>0 && extra>0)
                use=[use next];
                replace=replace-1;
                extra=extra-1;
                next=next+1;
            end;
            
                %cheack again that there is usabale coils information if not we won't fit this box 
            if length(use)>=4 
                useOLD=use;
                        %Fit the data this with  a nested function (below)
                [res(:,:,nn) resnorm(nn) exitflag(nn) num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt) ;
                
                        %check that the number of corect fitted coil are ther number that are in
                if length(use)==length(useOLD);
                    %all done we can move on; save the coil we end up fitting
                    
                    coilList(1:length(use),nn)=coils_list(use);
                    
                elseif num>=4,%
                    %last try, not using all coils but the one that seem
                    %fine
                                            %Fit the data this with  a nested function (below)
                    [res(:,:,nn) resnorm(nn) exitflag(nn) num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt);
                    if find(nopossible);
                        res(:,:,nn)=0;
                        skip(nn)=1;
                        disp('skipping this area - can not fit it')
                    else
                        %save the coil we end up fitting and move on
                        coilList(1:length(use),nn)=coils_list(use);
                    end;
                    
                else
                    %if every thing failed and 
                    %still there are too many coils need to be
                    %replace,  this is a bad box let move on with out
                    %fitting
                    res(:,:,nn)=0;
                    skip(nn)=1;
                    disp('skipping this area - can not fit it')
                end
            else
%                failed
                res(:,:,nn)=0;
                skip(nn)=1;
                disp('skipping this area - can not fit it')
            end
            %                 elseif num<5, %too many coils need to replace,  this is a bad box let move on
            %                     res(:,:,nn)=zeros(x0);
            %                     skip(nn)=1;
            %                     disp('skipping this area - can not fit it')
            
        end
        
        
    end;
    %next box
    nn=nn+1
    toc
    %%
end

%% Save
%the name of the saved file is usfule to detecet what box are in it (st ed)
name=[ opt{1}.name '_' num2str(st) '_' num2str(ed)];

save(name,'res','resnorm','exitflag','st','ed','skip' ,'coilList', 'str')
%saving the information about those boxes we fit.
%the coils fitted polynomial coefishents (res)
%the fittting residuals (resnorm)
%the exis flag that to end the fit  (exitflag)
%the box numbers (st ed)
%the boxes  we skip
% the coil list we used for we each boxs (coilList)
%Poly the polinomyal we fit the bais with



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%nested functions
%% I. Fit the box data

function [res resnorm exitflag num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt)

%initilaized the fit
%qoptions =  optimset('LevenbergMarquardt','on','Tolx',1e-6,'TolF',1e-6,'MaxIter',300,'MaxFunEvals',inf,'Display', 'iter');
a=version('-date');
if str2num(a(end-3:end))==2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'off','Tolx',1e-12);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'off','Tolx',1e-12);%'TolF',1e-12
    
end
inDat=inDat1(:,use);
bm=bm1(:,use);
x0=x01(use,:);
coefdat=tril(corrcoef(inDat),-1);


%call to lsq fit
[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par) errlocalGainUC_v2(par,inDat,opt{1}.Poly(bm(:,1),:),coefdat,bm,length(use)),x0,opt{1}.lb(1:length(use),:),opt{1}.ub(1:length(use),:),options);
res=zeros(opt{1}.numIn,size(res1,2));
res(1:size(res1,1),:)=res1;

%let do a qualety control by a secound  function. we check
%that the fit is fine  
[nopossible]= errlocalGainUC_vQualety(res,inDat,opt{1}.Poly,coefdat,bm,length(use),opt{1}.boxS);

%select the coil to use after check
use=use(find(nopossible==0));

num=length(use);
