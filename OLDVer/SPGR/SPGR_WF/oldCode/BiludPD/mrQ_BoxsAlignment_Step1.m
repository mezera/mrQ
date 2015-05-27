function [mat err matdc BOX bad]= mrQ_BoxsAlignment_Step1(M,opt,Fits,CoilsList,sqrtF)
%
% [mat err Reference matdc ]= mrQ_BoxsAlignment(M,opt,Fits,CoilsList,coils,sqrtF)
%the PD fit of each brain box need to be combined with the rest.
% all the boxes have a mean value of about 1
% that means that in order to combine the boxes (each might have a differet
% mean depend on where it is in the brain), we need to find a scalre to
% multipe each of the boxes.
% to find that will use the fact that there is overlap between the boxs so
% a scalre between each two overlap box can be calculate.
% when the box are scaled substracting two overlap region of two boxes
% should eqal to 1.
% we will use this idea
% we will build a matrix that estimate that have all the scalaer that used
% to and substract each two overlap box so it addes to zeros
%this matix can be soveld as alinear eqation system if one box set to be
%one.
%this box will be could the Reference box
% this function build a matrix mat that will alow us to estimate the scalre
% parameter of each box using this eqation
%  mat*C=y
% y is zeors (box X 1) vector with a single 1 in the reference box location
% mat is the eqtions of substract the niboring part of the boxs
% C is the scals of the box (that we will try to estimate globaly)
%
%we will loop over the boxes. and esine the values in the matrixes above in three step
%with three nested function  mrQ_getBoxsdat.
% 1. a functioon that load the box information and data and fits
% 2. we will find the feel the reference box values and it nighbors for
% the matrixes. mrQ_AvrageBox_reference
% 3. then we will add all the other boxes and nigbers. mrQ_getBoxsdat
%
% INPUTS:
%   opt         - a structure with the parameters that was used to optimaized the coil gains
%
%   M0          - The combined/aligned M0 data
%
%   Fits        - the coefitents that where fited to the polinomial that
%                 explain the coils gains (for eacjh box)
%
% CoilsList     - the list of coils that where use to fit (for each box)
%
%   sqrtF       - use the sqrt of the signal (zero mean no 1 yes)
%
% OUTPUTS
%
% mat            - the matrix we set to solve the scalers (see the
%                   description abouve)  (size (box number X box number) )
% err            - a matrix that estimate the miss overlap after scaling
%                   between each two boxs (size (box number X box number) )
%
% Reference      - the box we set to be ther refernce box (arbitratry any
%                   box with many useable nighboring boxs)
% matdc          - a matrix that of the scale term for each two boxs (size (box number X box number) )


%% CHECK INPUTS  SET DEFAULTS and intilaized parameters


list=opt{1}.wh;
skipf=0;
if(~exist('sqrtF','var')  || isempty(sqrtF))
    sqrtF=0;
end
SZ=size(opt{1}.X);
now=1;
done=0;
bad=[];
boxs=[];
reference=-1;
mat=zeros(length(list));
err=mat;
matdc=err;
mask= zeros(size(M.data(:,:,:,1)));

%
%% now lets repeat it for all the boxes it the list
done=zeros(length(list),1);
%go over the list of oxes
while now~=length(list)
    % get the location
    do_now=list(now);
    [fb(1,1) fb(1,2) fb(1,3)]=ind2sub(SZ,do_now);
    
    
    
    %1.let find and get the info from  all the  nigboring boxs as well
    boxs=[];
    boxlist=[];
    count=1;
    
    % loop over the near boxes
    for Xx= fb(1)-1:fb(1)+1
        for Yy= fb(2)-1:fb(2)+1
            for Zz= fb(3)-1:fb(3)+1
                % find the  box number
                
                if( Xx>SZ(1) || Yy>SZ(2) || Zz>SZ(3) || Xx<1 || Yy<1 || Zz<1)
                    %this location is not exsisting lets move on
                else
                    
                    boxNum=sub2ind(size(opt{1}.X),Xx,Yy,Zz);
                    %find where this box data is saved
                    Bwh=find(opt{1}.wh==boxNum);
                    
                    
                    if (~isempty(Bwh))
                        if (~any(bad==Bwh) )
                            if done(Bwh)==1
                                box=BOX{Bwh}.box;
                              
                                boxs{count}=box;
                                    boxs{count}.boxID=Bwh;
                                    boxlist(count)=Bwh;
                                    count=count+1;
                            else
                                %if this box is found to bad there is no point working on it. if it the Reference we don't what it as well.
                                fb1=[Xx Yy Zz];
                                % get the box fits information or skip it if it a bad
                                % fit
                                [skip box]=mrQ_getBoxsdat(fb1,opt{1},CoilsList(Bwh,:),Fits(:,:,Bwh),M,sqrtF,mask);
                                BOX{Bwh}.skip=skip;
                                if skip ~=1;
                                    %record the information
                                    boxs{count}=box;
                                    boxs{count}.boxID=Bwh;
                                    boxlist(count)=Bwh;
                                    count=count+1;
                                    BOX{Bwh}.box=box;
                                    
                                else
                                    % record the bad fits boxs
                                    bad=[bad Bwh];
                                end
                                done(Bwh)=1;
                            end
                        end
                    end
                end
            end
        end
    end
    % avaluate the matrixes for the  box we select and overlap  box around it
    [mat err matdc]=mrQ_AvrageBox(boxs,mat,err,boxlist,matdc);
    now=now+1;
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nested functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% function 1
%

function [skip box]=mrQ_getBoxsdat(fb,opt,coilIn,coilfit,M,sqrtF,mask)
%
%INPUTS
% 
% fb           - the box location in mashgrid of the the image space.
%
% opt           -  a structure with the parameters that was used to optimaized the coil gains
%
% coilIn        -  the list of coils that where use to fit the box
%
% coilfit        - the coefitents that where fited to the polinomial that
%                 explain the coils gains in the box
% 
%   M          - The combined/aligned M0 data
%%
%   mask       - an intilaiztion with zeros image that used to set the mask
%                   size(mask)=size(M(:,:,:,1))
%
%OUTPUTS
%
%boxs           - the structure of the nearby boxs with data mask and locations
%
% skip          - if somthing is wrong with this boxs data or fit skip = 1 and it won't
%                   be used 


%% CHECK INPUTS AND SET DEFAULTS

coils=(find(coilIn));
if isempty(find(coilIn))
    box=0;skip=1;
    return
end;
if isempty(find(coilfit))
    box=0;skip=1;
    return
end;
if ~isempty(find(isnan(coilfit)))
    box=0;skip=1;
    return
end;
%
%get the location
[Xx Yy Zz skip]=MrQPD_boxloc(opt,fb);
if skip==1;
    box=0;
    return
end;

%%
%load the raw M0 data that was fitted
if sqrtF==1
        %for the case of sqrt on theM0 images (not the defult)
    inDat=double(sqrt(M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),coilIn(coils))));
    
else
    inDat=double((M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),coilIn(coils))));
end
% get the fitted coefisent of the coil gain estimated by polynomials
Gain1=zeros(size(inDat));
Val=Gain1;
% get the box mask
bm=opt.brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));

%calculate PD (val1) from raw M0 images and coils gain
for i=1:size(inDat,4),%opt.numIn
    %
    Gain1(:,:,:,i) = reshape(opt.Poly*coilfit(i,:)',opt.boxS);
    Val(:,:,:,i) = inDat(:,:,:,i)./Gain1(:,:,:,i);
end;
%ResVal=median(Val,4);
% get the avrage PD fit of the different coils
ResVal=mean(Val,4); %

% box keeping of locations
box.Xx=Xx;box.Yy=Yy; box.Zz=Zz;

% mask= zeros(size(M.data(:,:,:,1)));
% mask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=1;
% mask(~logical(opt.brainMask))=0; %let work only when we fit  %%% take to long!!!

%make a mask location of the box in the imaging space
mask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=opt.brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
mask=logical(mask);

%save the location
box.wh=find(mask);
%box.wh=box.wh(:);




%%
%control for outlayers 
%lets wait by coils by signal
W=inDat; 
for i=1:size(inDat,4)
    t=inDat(:,:,:,i);
    W(:,:,:,i)=mean(t(:));
end
W=W./sum(W(1,1,1,:));

%ResVal=sum(Val.*W ,4); %waited the coils by SNR'


%get an estimation of the mean signal
MR=mean(ResVal(bm));

%stdE=std(ResVal,[],4);

%gert thebox  M0 data if we won't remoce the gain bias
Dat=sum(inDat.*W ,4);
Dat=(Dat./mean(Dat(bm))).*mean(ResVal(bm));
%stdI=std(Dat,[],4);
clear inDat


%ResVal=mean(Val,4);

%cheack for variation in pd (some time we get crazy high or low values this
%then std is high and we don't want the box or at list the few voxel in it
c=((std(Val,[],4)));
box.Val=ResVal(bm);
box.fb=fb;
% calculate the variance in the box before and after coils gain removed. in
% the normal case removing the gain will decrise the variationin the image
box.SD=std(ResVal(bm));
box.SDI=std(Dat(bm));

SS=1-(ResVal(bm)-c(bm))./ResVal(bm);
%SS=c(bm);
%so let check for crazy fit and exclude them
if (box.SD<box.SDI.*1.1); %we first check that this is aresnable fit if the std of the box grow this is probabaly becouse insate of a fit the values expload in that case we won't use this box
    
    if (any(box.Val<0)) %if we still have few nagative we won't use them for alighnment (happan in the edge of the brain air or noise voxels
        wh=find(box.Val>0);
        box.Val=box.Val(wh);
        box.wh=box.wh(wh);
        SS=SS(wh);
    end
    if any(box.Val>(MR+3*box.SD)) %if we still have few very high values we won't use them for alighnment (happan in the edge of the brain air or noise voxels or some csf voxel that have very low SNR)
        wh=find(box.Val<(MR+3*box.SD));
        box.Val=box.Val(wh);
        box.wh=box.wh(wh);
        SS=SS(wh);
        
    end
    if  any(box.Val<(MR-3*box.SD))%if we still have few very low values we won't use them for alighnment (happan in the edge of the brain air or noise voxels or some csf voxel that have very low SNR)
        wh=find(box.Val>(MR-3*box.SD));
        box.Val=box.Val(wh);
        box.wh=box.wh(wh);
        SS=SS(wh);
        
    end
    
    if  any(SS>0.06)% if part of this box is unconclusive so the std between the different coils is very high we better not use it. that happean it the edge of the boxs becouse of miss fit; or when single to noise is low  like csf or air edge
        wh=find(SS<0.06);
        box.Val=box.Val(wh);
        box.wh=box.wh(wh);
    end
    box.Median=median(box.Val(:));
     box.Val=box.Val./median(box.Val(:));
        
    
else
    skip =1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%% function 2
%

function [mat err matdc]=mrQ_AvrageBox(boxs,mat,err,boxlist,matdc)
%we will find the scale of each box so it will be the same with the values in the
%overlap boxes ( a scalar term). the box are free from gain and the only different is a scaler that need to be fit.
%we build the matrix mat that sumerized all of those scaler for all the boxes. the bix mat  will alow us to find the align values for the
%all boxes thogther by a minization of all the scalares thogther. the
%values in mat are added thogter so the fitting problem can be solved. the
%idea the box(i) -box(j)*scaler equal 0. we will build this eqation for all
%the overlap boxes and solve them thogther in the end. (out side of this
%function)
%we will cheack that the align error after alignment with the scaler.
%(big error mean bad aligment and therefor bad fit of one of the boxes). this is saved in an err
%matrix that is in the same size of mat. the term err(i,j) is the error
%between box i and j after the scaler alignment
%output the mat and err of after adding the term of the current boxes.
%matdc is the scale given for each to boxes and it is the same size as mat
%INPUTS
% boxs           - the structure of the nearby boxs with data mask and locations
%
% mat            - the matrix we set to solve the scalers (see the
%                   description abouve)  (size (box number X box number) )
%                   the new boxes information added to it
% err            - a matrix that estimate the miss overlap after scaling
%                   between each two boxs (size (box number X box number) )
%                   the new boxes information added to it
%
% boxlist        - the list of boxes we are using now
% Reference-     - the box that was set to be refernce box (we will avoide don
%                     to use it here)
%
% matdc          - a matrix that of the scale term for each two boxs (size (box number X box number) )
%                   the new boxes information added to it
%
%OUTPUTS
%
% mat            - the matrix we set to solve the scalers (see the
%                   description abouve)  (size (box number X box number) )
% err            - a matrix that estimate the miss overlap after scaling
%                   between each two boxs (size (box number X box number) )
%
% matdc          - a matrix that of the scale term for each two boxs (size (box number X box number) )


Buse=1:length(boxs);


for jj=Buse %run over the boxes
    
    
        BBuse=Buse(Buse~=jj);
        for j=1:length(BBuse)
            i=BBuse(j);% run over  a box vrs. the other boxes
            
            %if  mat(i,jj)~=1           %if this boxes was already calculated before we won't do it again
            if  mat(boxs{i}.boxID,boxs{jj}.boxID)~=1           %if this boxes was already calculated before we won't do it again
                
                [tf, loc]=ismember(boxs{i}.wh,boxs{jj}.wh); % the overlap locations
                loc1=loc(loc>0);
                
                if ( length(loc1)/length(boxs{jj}.wh)>0.1 && length(loc1)/length(boxs{i}.wh)>0.1)  % check if there is at list 10% overlap to calculate the scaler
                    
                    dc=median(boxs{jj}.Val(loc1)./(boxs{i}.Val(tf)));  %the scaler
                    
                    % err(i,jj)=median(abs(boxs{jj}.Val(loc1)-(boxs{i}.Val(tf)).*dc)./((boxs{jj}.Val(loc1)+(boxs{i}.Val(tf)).*dc)./2)); %the error
                    err(boxs{i}.boxID,boxs{jj}.boxID)=median(abs(boxs{jj}.Val(loc1)-(boxs{i}.Val(tf)).*dc)./((boxs{jj}.Val(loc1)+(boxs{i}.Val(tf)).*dc)./2)); %the error
                    if (err(boxs{i}.boxID,boxs{jj}.boxID)<0.05 && dc>0.5 && dc<2) %this is too much they don't fit realy

                    %if (err(boxs{i}.boxID,boxs{jj}.boxID)<0.05) %this is too much they don't fit realy
                        % update of the scaler matrix
                        
                        mat(boxs{i}.boxID,boxs{jj}.boxID)=1;
                        mat(boxs{i}.boxID,boxs{i}.boxID)=mat(boxs{i}.boxID,boxs{i}.boxID)+(-dc);
                        matdc(boxs{i}.boxID,boxs{jj}.boxID)=dc;
                        
                        %                         if dc>1.5 || dc<0.5
                        %                             keyboard
                        %                         end
                        %                 mat(i,jj)=1;
                        %                 mat(i,i)=mat(i,i)+(-dc);
                    end
                else
                    %keyboard
                end
            end
        end
    
end

