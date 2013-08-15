function [M0,Wait,VxUsed,Boxorder,err1]= mrQ_BoxsStepWiseAlignment(BOX,matdc,bad,boxNm,opt)


DoNow=boxNm;
%DoNext=boxNm;
SZ=size(opt{1}.X);
ToDo=ones(size(opt{1}.wh));
ToDo(bad)=0;
TodoNext=zeros(size(opt{1}.wh));
doneSt=0;
M0=zeros(size(opt{1}.brainMask));
Wait=M0;
count=1;
while any (ToDo==1)
tic
%for i=1:length(DoNext)


wh=find(matdc(:,DoNow)); wh=wh';
wh=[DoNow wh];


 MMi=nan(length(opt{1}.brainMask(:)),length(wh));
 Mdc=matdc(wh,DoNow);
Mdc(1)=1;
 for i=1:length(wh)
     MMi(BOX{wh(i)}.box.wh(:),i)=BOX{wh(i)}.box.Val.*Mdc(i);
 end
 Im1=sum(~isnan(MMi),2);
 use=find(Im1>6);
 
MMii=nan(length(use),length(wh));
MMii(:,:)= MMi(use,:);
clear MMi
ImS=nanstd(MMii,[],2);
 Im=nanmedian(MMii,2);
 good=(ImS./Im)<0.03;
 Im=Im(good);
 use=use(good);
  
 VxUsed(DoNow,1)=length(use);
  
 if doneSt==0
     doneSt=1;
     M0(use)=Im;
     Wait(use)=1;
     VxUsed(DoNow,2)=VxUsed(DoNow,1);
     err1(DoNow)=0;
 else
     whV=find(M0);
     
     [tf, loc]=ismember(use,whV); % the overlap locations
 loc=loc(loc>0);
 if ~isempty(loc)
 dc=median(M0(whV(loc))./Im(tf));
 VxUsed(DoNow,2)=length(find(tf==0));
      err1(DoNow)=median( (M0(whV(loc))-Im(tf).*dc)./M0(whV(loc)));

     M0(use)=(M0(use).*Wait(use)+Im.*dc)./(Wait(use)+1); %waited sum
     Wait(use)=Wait(use)+1;
 end
 end
 %book keeping
 Boxorder(DoNow)=count;
 count=count+1
 


 
ToDo(DoNow)=0;

TodoNext(wh(2:end))=1;
TodoNext(ToDo==0)=0;
clear Im1 MMii use Im good ImS wh Mdc tf loc dc whV
% lets find next by maxsimum overlap
wh=find(TodoNext);
whV=find(M0);
for i=1:length(wh)
     tf =ismember(BOX{wh(i)}.box.wh(:),whV); % the overlap locations
               overlap(i)=length(find(tf));
end
if notDefined('overlap')
    toTry=find( ToDo==1);
    if isempty(toTry)
   ToDo(:)=0;
    else
     DoNow(toTry(1));   
    end
else
    mostOverLap=find(overlap==max(overlap));
    DoNow=wh(mostOverLap(1));
    clear wh whV tf  overlap
    toc
end

end




