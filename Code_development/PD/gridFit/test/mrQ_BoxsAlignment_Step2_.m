function [mat,A]= mrQ_BoxsAlignment_Step2(BOX,err,matdc,bad,rounds,boxNm,opt)


SZ=size(opt{1}.X);
B=opt{1}.wh(boxNm);
[fb(1,1) fb(1,2) fb(1,3)]=ind2sub(SZ,B);

BB=boxNm;

if any(bad==BB)
    return
    A=[];
    mat=[];
end
%1.let find and get the info from  all the  nigboring boxs as well
boxs=[];
boxlist=[];
count=1;
minWH=prod(size(opt{1}.brainMask));
maxWH=0;

wh=find(matdc(:,BB)); wh=wh';
wh=[BB wh];
done=0;
wh1=wh;
round=0;
while done==0

 for i=1:length(wh); wh1=[wh1 find(matdc(:,wh(i)))']; end;
wh=unique(wh1);
wh1=wh;
round=round+1
 if round==rounds
     done=1;
 end

end;

for i=1:length(wh)
     Bwh=wh(i)
     boxs{count}=BOX{wh(i)}.box;
     boxs{count}.boxID=Bwh;
     boxlist(count)=Bwh;
     minWH=min(minWH,min(boxs{count}.wh));
     maxWH=max(maxWH,max(boxs{count}.wh));
     count=count+1;
end



 
  reference=find(boxlist==boxNm);

[mat]=mrQ_AvrageBox_reference(boxs,reference,matdc,err);




y=zeros(size(mat,1),1);
y(reference)=1;
%solve it as multi  linear eqation
C=pinv(mat'*mat)*mat'*y;
% the C we need is one over the one we fit
C1=1./C;

st=1;
jump=minWH-1;
ed=maxWH-jump;
A=nan(ed,count-1);
wh2=find(C1>0.1 & C1<2);

 B=zeros([size(opt{1}.brainMask) 2]);
 C=zeros(size(opt{1}.brainMask));
 Bnm=[55 54];
 for i=1:2

 C=zeros(size(opt{1}.brainMask));

 C(boxs{Bnm(i)}.wh)=boxs{Bnm(i)}.Val;
 B(:,:,:,i)=C;
 end
% 
% 
% for i=1:length(wh2)
%     cc =wh2(i);
%     A(boxs{cc}.wh(:)-jump,i)=boxs{cc}.Val.*C1(i);
% end
% a=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%% function 1
%
function [mat]=mrQ_AvrageBox_reference(boxs,reference,matdc,err)
%we will find the scale of each box so it will be the same with the values in the
%overlap boxes ( a scalar term). the box are free from gain and the only different is a scaler that need to be fit.
%we build the matrix mat that sumerized all of those scaler for all the boxes. the bix mat  will alow us to find the align values for the
%all boxes thogther by a minization of all the scalares thogther. the
%values in mat are added thogter so the fitting problem can be solved. the
%idea the box(i) -box(j)*scaler equal 0. we will build this eqation for all
%the overlap boxes and solve them thogther in the end.(out side of this
%function)
% in the case of the reference we deside that  mean(box(refence))*scaler=1;
%the refernce box is unique as instade of adding it to the other box we say
%it eqal an arbitrariy value (1). so all the boxes will be align to its mean values.
%we will cheack that the align error after alignment with the scaler.
%(big error mean bad aligment and therefor bad fit of one of the boxes). this is saved in an err
%matrix that is in the same size of mat. the term err(i,j) is the error
%between box i and j after the scaler alignment
%output the mat and err of after adding the term of the current boxes.
%matdc is the scale given for each to boxes and it is the same size as mat
%INPUTS
% boxs           - the structure of the nearby boxs with data mask and locations
%
% reference-     - the box that was set to be refernce box
%
% mat            - the matrix we set to solve the scalers (see the
%                   description abouve)  (size (box number X box number) )
%                   the new boxes information added to it
% err            - a matrix that estimate the miss overlap after scaling
%                   between each two boxs (size (box number X box number) )
%                   the new boxes information added to it
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
%Buse=Buse(Buse~=reference);
mat=zeros(length(boxs));
mat(reference,reference)=1;

for jj=Buse % run over the boxes
    BBuse=Buse(Buse~=jj);
    for j=1:length(BBuse) % run over  a box vrs. the other boxes
        i=BBuse(j);
        if i~=reference
        if matdc(boxs{i}.boxID,boxs{jj}.boxID)~=0
            dc=matdc(boxs{i}.boxID,boxs{jj}.boxID);
            
             mat(i,jj)=1;
            
            if jj==reference                
                mat(i,i)=-dc;
            else
                
                mat(i,i)=mat(i,i)+(-dc);
            end
            
        end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


