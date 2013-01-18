function mrQfitPD_multiCoils_M0(outDir,M0cfile,degrees)


if(~exist('degrees','var'))
    disp(['we will use the defult polynomyals degree is of 2 for coil estimation']);
    degrees=2;
end

if(~exist('SunGrid','var'))
    disp(['no SunGrid use']);
    SunGrid=0;
end

if(~exist('M0cfile','var'))
    M0cfile{1}=[outDir '/AligncombineCoilsM0.nii.gz'];
end


logname=[outDir '/fitLog.mat'];
 load(logname);

%opt{1}.degrees=degrees;
%opt{1}.outDir=outDir;
%numIn=10;


BMfile=fullfile(outDir,'brainMask.nii.gz');

if(exist(BMfile,'file'))
    disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    xform=brainMask.qto_xyz;
    mmPerVox=  brainMask.pixdim;
    brainMask=logical(brainMask.data);

else
    disp(['error , can not find the file: '   BMfile]);
end;

cdffileFS=fullfile(outDir,['csf_FS_T1.nii.gz']);
if(exist(cdffileFS,'file'))
    M0f = readFileNifti(cdffileFS);
    M0f=double(M0f.data);
    M0f=M0f.*100;
else
    disp(['error , can not find the file: '   cdffileFS]);
    error
end;

[a b]=fileparts(fileparts(fileparts(fileparts(fileparts(outDir)))));
sb=b;
disp([' subject name for M0 fit is ' b]);

Avmap=zeros(size(brainMask));
M0ft=zeros(size(brainMask));
M0ftfull=zeros([size(M0ft) 27]);
%% the mask we will use to find the boxes to fit
NHOOD=[0 1 0; 1 1 1; 0 1 0];
NHOOD1(:,:,2)=NHOOD;

NHOOD1(:,:,1)=[0 0 0; 0 1 0; 0 0 0];
NHOOD1(:,:,3)=NHOOD1(:,:,1);
SE = strel('arbitrary', NHOOD1);

%bb=imdilate(b,SE);
%%
%we will try to till the brain with boxes of rafly 20mm^3 and with overlap
%of 2;
% sz=(size(brainMask));
% boxS=round(24./mmPerVox);
% even=find(mod(boxS,2)==0);
% 
% boxS(even)=boxS(even)+1;
% opt{1}.boxS=boxS;
% %overlap=round(boxS.*0.1); %10% over lab
% %overlap=round(boxS.*0.7); %70% over lap
% overlap=round(boxS.*0.5); %70% over lap
%the center of the CSF will be the first box to fit

% the center of the boxs we will fit
%[opt{1}.X,opt{1}.Y,opt{1}.Z] = meshgrid(round(boxS(1)./2):boxS(1)-overlap(1):sz(1)    ,    round(boxS(2)./2):boxS(2)-overlap(2):sz(2)    , round(boxS(3)./2):boxS(3)-overlap(3):sz(3));


controlmask=zeros(size(opt{1}.X));
%donemask=controlmask;
donemask=opt{1}.donemask;
%opt{1}.HboxS=(boxS-1)/2;

%two mask that control what box done and which need to e done





%box=CSF(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3));




%%


%let first find on what box we will work
%  tmpfile=fullfile(outDir,'tmpM0fit.mat');
%     if (exist(tmpfile,'file'))
%         load (tmpfile)
%     else


for i=1:length(opt{1}.wh),%run over the box you like to fit and cheack there is data there

    [fb(1) fb(2) fb(3)]=ind2sub(size(opt{1}.X),opt{1}.wh(i));

    [csf(i)]= mrQM0fit_isBoxcsf(opt{1},brainMask,fb,M0f) ;

    
end;
[xs ys zs]=ind2sub(size(opt{1}.X),opt{1}.wh(find(csf==max(csf))));

controlmask(xs(1),ys(1),zs(1))=1;



jumpindex=opt{1}.jumpindex;
%jumpindex=10;
jobindexs=1:ceil(length(opt{1}.wh)/jumpindex);

for j=1:length(M0cfile)
    Fits=zeros(opt{1}.numIn,jumpindex,length(opt{1}.wh));
    resnorms=zeros(length(opt{1}.wh),1);
    exitflags=zeros(length(opt{1}.wh),1);
    CoilsList=zeros(length(opt{1}.wh),opt{1}.numIn);
    %dirname=[outDir '/tmpSGM0' num2str(j)];

    opt{1}.brainMask=brainMask;

    JBname=opt{1}.name;%[dirname '/M0boxfit_iter'] ;
    %%

    for i=1:length(jobindexs);



        st=1 +(jobindexs(i)-1)*jumpindex;
        ed=st+jumpindex-1;
        if ed>length(opt{1}.wh), ed=length(opt{1}.wh);end;

        Cname=[JBname '_' num2str(st) '_' num2str(ed)];
        load(Cname);
        ss=length(resnorm);

        Fits(:,:,st:st+ss-1)=res;
        resnorms(st:st+ss-1)=resnorm;
        exitflags(st:st+ss-1)=exitflag;
        
        CoilsList(st:st+ss-1,:)= coilList(:,1:ss)';

    end;
%%

    M=readFileNifti(M0cfile{j});
    coils=size(M.data,4);
%profile on
    done=0;tt=1;
    while done==0

        do_now=find(controlmask==1);
        if isempty(do_now)
            done=1;
        else

            for kk=1:length(do_now)
                BB=find(opt{1}.wh==do_now(kk));
                if isempty(BB)
                                    keyboard;

                 %   ha=1;
                else
                end
%     if tt==155
%                                      keyboard;
%                  end
    
                [fb(1,1) fb(1,2) fb(1,3)]=ind2sub(size(controlmask),do_now(kk));
                if tt==1,
                    [known,inDat,skip use ]= mrQM0_Partofit(opt{1},M,fb,coils,M0f,CoilsList(BB,:));
                else
                    [known,inDat,skip use]= mrQM0_Partofit(opt{1},M,fb,coils,M0ft,CoilsList(BB,:));
                end;
                tt=tt+1
                G=Fits(use,:,BB);
%                    if find(CoilsList(BB,:)==0);
%                                                               keyboard;
%                    end
% %                   tt=tt+1
%                  if tt==848
%                                      keyboard;
%                  end
                if (skip==0  && isempty(find(isnan(G))) && ~isempty(find(G)) && ~isempty(find(CoilsList(BB,use))) ) ;
%                     if exitflags(BB)~=0
%                     keyboard;
%                     end
%                      if resnorms(BB)<10
%                     keyboard;
%                      end
%                     if resnorms(BB)>3000
%                  
%                      end
                [M0f,donemask,Avmap,M0ft,M0ftfull] =mrQM0fiti_Res(opt{1},fb,G,inDat,M0f,donemask,known,Avmap,M0ft,M0ftfull,resnorms(BB),exitflags(BB));
                else
                end;
            end
        end

        controlmask=imdilate(controlmask,SE);
        controlmask=controlmask+donemask;


    end;
    
    
    
    
         WFfile=fullfile(outDir,['PDsqrt_fitbox_' num2str(j) '.nii.gz']);

 dtiWriteNiftiWrapper(single(M0f), xform, WFfile);
%                keyboard;
PDmedianfromGain(M0ftfull,xform,outDir,j);


    end

%%%%%%%%%%%%%%%
