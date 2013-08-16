function[known,inDat,skip coils]= mrQM0_Partofit(opt,M,fb,coils,M0,In)

coils=(find(In));
if isempty(find(In))
   known=0;inDat=0;skip=1;
    return
end; 


[Xx Yy Zz,skip]=MrQPD_boxloc(opt,fb);

if skip==1;
   known=0;inDat=0;
    return
end;



box(:,:,:,:)=double(M.data(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2),:));


% for i=1:coils;rate(i)=mean(mean(mean(box(:,:,:,i))));end;
% [v In]=sort(rate,'descend');

tt=M0(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
bm=opt.brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));

for i=1:length(coils),
    tmp=box(:,:,:,In(i));
    no=find(isnan(tmp));
    bm(no)=0;
    no=find(isinf(tmp));
    bm(no)=0;
    Mm=median(tmp(bm));
    Sd=std(tmp(bm));

    no=find(tmp>(Mm+3*Sd));
    bm(no)=0;
    no=find(tmp<(Mm-3*Sd));
    bm(no)=0;
    no=find(tmp<=0);
    bm(no)=0;


end;

if length(find(bm))<20;
   known=0;inDat=0;skip=1;
    return
end;

tt(~bm)=0;

if find(tt)
known(1,:)=find(tt);
known(2,:)=tt(known(1,:));
else
 known(1,2)=0;   
end

inDat=double(sqrt(box(:,:,:,In(coils))));
%inDat=double(sqrt(box(:,:,:,In(1:opt.numIn))));