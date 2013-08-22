    inDir='/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/SPGR_2/Align_0.9375_0.9375_1'
M0cfile = fullfile(inDir,'AligncombineCoilsM0.nii.gz');


 load('/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/test/fitLog_corect1.mat');
outDir1='/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/test/tmpSGM01_corect1';

opt{1}.name='/biac4/wandell/biac2/wandell2/data/WH/008_AM/Qmr/20111020_1294_32ch_1mm3/20111020_1294/test/tmpSGM01_corect1/M0boxfit_iter'

JBname=opt{1}.name;

if isfield(opt{1},'sqrtF')
    sqrtF=(opt{1}.sqrtF);
end
if(~exist('sqrtF','var')  || isempty(sqrtF))
    sqrtF=0;
end


jumpindex=opt{1}.jumpindex;
jobindexs=1:ceil(length(opt{1}.wh)/jumpindex);

BMfile=fullfile(inDir,'brainMask.nii.gz');

if(exist(BMfile,'file'))
    %  disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    xform=brainMask.qto_xyz;
    mmPerVox=  brainMask.pixdim;
    brainMask=logical(brainMask.data);
    
else
    disp(['error , can not find the file: '   BMfile]);
end;


%%


%% load the fit

[opt{1}.Poly,opt{1}.str] = constructpolynomialmatrix3d(opt{1}.boxS,find(ones(opt{1}.boxS)),opt{1}.degrees);

%initilaized the parameters
Fits=zeros(opt{1}.numIn,size(opt{1}.Poly,2),length(opt{1}.wh));
resnorms=zeros(length(opt{1}.wh),1);
exitflags=zeros(length(opt{1}.wh),1);
CoilsList=zeros(length(opt{1}.wh),opt{1}.numIn);

opt{1}.brainMask=brainMask;
%for ST=1:length(opt{1}.wh);strs{ST}=[];end

% loop over the fitted box files and load the parameters to  matrixs
for i=1:length(jobindexs);
    
    
    
    st=1 +(jobindexs(i)-1)*jumpindex;
    ed=st+jumpindex-1;
    if ed>length(opt{1}.wh), ed=length(opt{1}.wh);end;
    
    Cname=[JBname '_' num2str(st) '_' num2str(ed)];
    load(Cname);
    ss=length(resnorm);
    if isempty(res)
        res=0;
    end
    if ~strcmp(opt{1}.str,str)
        Cname
        error
    end
    Fits(:,:,st:st+ss-1)=res;
    resnorms(st:st+ss-1)=resnorm;
    exitflags(st:st+ss-1)=exitflag;
    
    CoilsList(st:st+ss-1,:)= coilList(:,1:ss)';
    
end;






%%
% load the M0 imaged that was used to the fit (multi coils 4D)
M=readFileNifti(M0cfile);
coils=size(M.data,4);

tic;[mat err matdc BOX bad]= mrQ_BoxsAlignment_Step1(M,opt,Fits,CoilsList,sqrtF);toc;

name=[outDir1 '/Step1_Fit' date];
    %toc
    save(name,'mat', 'err', 'matdc', 'BOX', 'bad')

[matBOX,A]= mrQ_BoxsAlignment_Step2(BOX,err,matdc,bad,1,773,opt);




[M0,Wait,VxUsed,Boxorder,err1]= mrQ_BoxsStepWiseAlignment(BOX,matdc,bad,773,opt);

name=[outDir1 '/Step2_Fit' date];
save(name,'M0','Wait','VxUsed','Boxorder','err1')





return


%%


%  [mat err Reference matdc]= mrQ_BoxsAlignment(M,opt,Fits,CoilsList,sqrtF);
%     
%     name=[outDir1 '/tmpFit_1' date];
%     %toc
%     save(name,'mat', 'err', 'Reference', 'matdc')
%     
%     
%     
%     
%     %%
%     %y=mat(Reference,:);
% %y=y';
% y=zeros(size(err,1),1);
% y(Reference)=1;
% %solve it as multi  linear eqation
% C=pinv(mat'*mat)*mat'*y;
% % the C we need is one over the one we fit
% C1=1./C;
% 
% % now when we have the scales we can combine the boxes
% %we like to exlude crazy boxs (very high or low scale
% % we know that C sepuuse to be around 1 if every thing went right
% 
% wh1=find(C1>0.1 & C1<2);
% %
% %in orther to avrage we need to keep track of the nmber of estimation of
% %each voxel (overlaping boxes makes it grater then one). we will keep
% %records on that by Avmap
% Avmap=zeros(size(M.data(:,:,:,1)));
% %we will save the avrage in M0f
% M0f=zeros(size(M.data(:,:,:,1)));
% %we will save the values for median avrage in M0full
% M0full=M0f(:);
%     
%     %%
%     for jj=1:length(wh1),
%     % jj
%     clear BB do_now fb Xx Yy Zz skip inDat In use  G Gain1 Val W t ResVal wh whh mask bm c SS  SD MR Val1
%     %what box to use
%     BB=(wh1(jj));
%     do_now=opt{1}.wh(BB);
%     
%     %set mask
%     mask=zeros(size(M.data(:,:,:,1)));
%     %get the location of the box in x,y,z cordinate
%     [fb(1,1) fb(1,2) fb(1,3)]=ind2sub(size(opt{1}.X),do_now);
%     [Xx Yy Zz,skip]=MrQPD_boxloc(opt{1},fb);
%     
%     %get the coil list we used in the box
%     In=CoilsList(BB,:);
%     use=(find(In));
%     use=(find(CoilsList(BB,:)));
%     
%     %load the raw M0 data that was fitted
%     if sqrtF==1
%         %for the case of sqrt on theM0 images (not the defult)
%         inDat(:,:,:,:)=double(sqrt(M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),In(use))));
%     else
%         inDat(:,:,:,:)=double(M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),In(use)));
%     end
%     % get the fitted coefisent of the coil gain estimated by polynomials
%     G=Fits(use,:,BB);
%     
%     %calculate PD (val1) from raw M0 images and coils gain
%     for i=1:size(inDat,4),%opt.numIn
%         %
%         Gain1(:,:,:,i) = reshape(opt{1}.Poly*G(i,:)',opt{1}.boxS);
%         Val1(:,:,:,i) = inDat(:,:,:,i)./Gain1(:,:,:,i);
%     end;
%     
%     % we can wait the PD by SNR% we desice not to do that becouse it can bias
%     % the fits
%     % W=inDat; %lets wait by coils
%     % for i=1:size(inDat,4)
%     %     t=inDat(:,:,:,i);
%     %     W(:,:,:,i)=mean(t(:));
%     % end
%     % W=W./sum(W(1,1,1,:));
%     %ResVal=sum(Val1.*W ,4); %waited the coils by SNR
%     
%     % get the avrage PD fit of the different coils
%     ResVal=mean(Val1,4); %
%     
%     
%     % get the brain mask of the boxs in box space
%     bm=opt{1}.brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
%     % get the brain mask of the boxs in full imaging space
%     mask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=opt{1}.brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
%     mask=logical(mask);
%     
%     %control for outlayers
%     c=((std(Val1,[],4)));
%     wh=find(mask);
%     SS=c(bm);
%     Val=ResVal(bm);
%     SD=std(ResVal(bm));
%     MR=mean(ResVal(bm));
%     
%     if (any(Val<0)) %if we still have few nagative values we won't use them for alighnment (happan in the edge of the brain air or noise voxels
%         whh=find(Val>0);
%         Val=Val(whh);
%         wh=wh(whh);
%         SS=SS(whh);
%     end
%     if any(Val>(MR+3*SD)) %if we have very high values e won't use them (happan in the edge of the brain air or noise voxels or some csf voxel that have very low SNR)
%         whh=find(Val<(MR+3*SD));
%         Val=Val(whh);
%         wh=wh(whh);
%         SS=SS(whh);
%         
%     end
%     if  any(Val<(MR-3*SD))%if we still have few very low value (happan in the edge of the brain air or noise voxels or some csf voxel that have very low SNR)
%         whh=find(Val>(MR-3*SD));
%         Val=Val(whh);
%         wh=wh(whh);
%         SS=SS(whh);
%         
%     end
%     
%     if  any(SS>0.06)% if part of this box is unconclusive so the std between the different coils is very high we better not use it. that happean it the edge of the boxs becouse of miss fit; or when single to noise is low  like csf or air edge
%         whh=find(SS<0.06);
%         Val=Val(whh);
%         wh=wh(whh);
%     end
%     
%     %add this box data to the other
%     
%     % mutipal the result by it scaler
%     ResVal=ResVal.*C1(BB);
%     % for mean avraging
%     M0f(wh)=M0f(wh)+Val.*C1(BB);
%     Avmap(wh)=Avmap(wh)+1;
%     
%     %this is to mesure the median avraging
%     szz1=size(M0full,2);
%     Raw=max(max(Avmap(:)),szz1);
%     Col=length(M0f(:));
%     szz=[Col,Raw ];
%     
%     tmp=zeros(szz);
%     tmp(:,1:szz1)=M0full;
%     
%     wh0=  sub2ind(szz,wh,Avmap(wh));
%     tmp(wh0)=Val.*C1(BB);
%     M0full=tmp;
%     
% end
% %%
% % mean avrage the PD values
% M0f(find(M0f))=M0f(find(M0f))./Avmap(find(M0f));
% 
% %median avrage the PD values
% M0full(M0full==0)=nan;
% M0full=nanmedian(M0full,2);
% M0full=reshape(M0full,size(M0f));
% 
% 
% 
% PDfile2=fullfile(outDir1,['PD_fitGboxMedian_1.nii.gz']);
% 
% PDfile1=fullfile(outDir1,['PD_fitGboxmean_1.nii.gz']);
% dtiWriteNiftiWrapper(single(M0f), xform, PDfile1);
%     dtiWriteNiftiWrapper(single(M0full), xform, PDfile2);
%     
%     