function mrQ_FitM0_SGE_v1(opt,jumpindex,jobindex)
% 
% FitM0_sanGrid_v1(opt,jumpindex,jobindex)
% 
% Comments
% 
% INPUTS:
%       opt - 
%       jumpindex - 
%       jobindex  -
% 
%      
% 
% (C) Stanford University, VISTA
% 
% 

%%

maxnumcompthreads(1)
j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;
if ed>length(opt{1}.wh), ed=length(opt{1}.wh);end;

M=readFileNifti(opt{1}.dat);
coils=size(M.data,4);
Mall=mean(M.data,4);
%this is away to clear out layer but it can make problem so we can think
%about it again ...
Mall(Mall<mean(Mall(opt{1}.brainMask))-3*std((Mall(opt{1}.brainMask))))=0;
Mall=logical(Mall);
for i=1:coils,
    tmp=M.data(:,:,:,i);
    tmp(~Mall)=0;
    M.data(:,:,:,i)=tmp;
end;
%options =  optimset('LevenbergMarquardt','on','Tolx',1e-6,'TolF',1e-6,'MaxIter',300,'MaxFunEvals',inf,'Display', 'iter');
nn=1;

coilList=zeros(opt{1}.numIn,ed-st+1);
for i= st:ed,

    %run over the box you like to fit

    tic
    [fb(i,1,1) fb(i,1,2) fb(i,1,3)]=ind2sub(size(opt{1}.X),opt{1}.wh(i));

    %[x0,known,inDat,inNan]= mrQM0fitini_Param(opt{1}.X,opt{1}.Y,opt{1}.Z,M,fb(i,:),opt{1}.HboxS,opt{1}.boxS,coils,opt{1}.M0f,opt{1}.numIn,opt{1}.Poly,opt{1}.degrees,[],[],opt{1}.outDir,[],opt{1}.donemask,opt{1}.brainMask);
    [x01,known,inDat1,bm1,szB,skip(nn) coils_list]= mrQM0fitini_Param_v(opt{1},fb(i,:),M,coils);


    %use should define the brain mask --> and we should get outlayer out.
    if skip(nn)==0
        use=[1:min(size(inDat1,2),opt{1}.numIn)];

        [res(:,:,nn) resnorm(nn) exitflag(nn) num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt);
        if num==opt{1}.numIn
            %all done we move on
            coilList(1:length(use),nn)=coils_list(use);
        else%if num<5,
            %some coils are better replaced

            replace=opt{1}.numIn-num;
            extra=size(inDat1,2)-opt{1}.numIn;
            next=opt{1}.numIn+1;
            while (replace>0 && extra>0)
                use=[use next];
                replace=replace-1;
                extra=extra-1;
                next=next+1;
            end;

            if length(use)>=4
                useOLD=use;
                [res(:,:,nn) resnorm(nn) exitflag(nn) num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt) ;


                if length(use)==length(useOLD);
                    %all done we move on
                    coilList(1:length(use),nn)=coils_list(use);

                elseif num>=4,%
                    %last try, not using all coils

                    [res(:,:,nn) resnorm(nn) exitflag(nn) num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt);
                    if find(nopossible);
                        res(:,:,nn)=0;
                        skip(nn)=1;
                        disp('skipping this area - can not fit it')
                    else
                        coilList(1:length(use),nn)=coils_list(use);
                    end;

                else
                    %still too many coils need to
                    %replace,  this is a bad box let move on
                    res(:,:,nn)=0;
                    skip(nn)=1;
                    disp('skipping this area - can not fit it')
                end
            else
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
    nn=nn+1
    toc

end








name=[ opt{1}.name '_' num2str(st) '_' num2str(ed)];
save(name,'res','resnorm','exitflag','st','ed','skip' ,'coilList')












function [res resnorm exitflag num use nopossible] =FitDat(nn,inDat1,bm1,x01,use,opt)

options =  optimset('LevenbergMarquardt','on','Tolx',1e-6,'TolF',1e-6,'MaxIter',300,'MaxFunEvals',inf,'Display', 'iter');
inDat=inDat1(:,use);
bm=bm1(:,use);
x0=x01(use,:);
coefdat=tril(corrcoef(inDat),-1);


[res1, resnorm,dd1,exitflag] = lsqnonlin(@(par) errlocalGainUC_v1(par,inDat,opt{1}.Poly,coefdat,bm,length(use)),x0,opt{1}.lb,opt{1}.ub,options);
res=zeros(opt{1}.numIn,size(res1,2));
res(1:size(res1,1),:)=res1;
%let do a qualety control:

[nopossible]= errlocalGainUC_vQualety(res,inDat,opt{1}.Poly,coefdat,bm,length(use),opt{1}.boxS);

use=use(find(nopossible==0));

num=length(use);
