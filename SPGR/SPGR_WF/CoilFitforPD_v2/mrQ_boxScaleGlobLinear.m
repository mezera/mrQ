function [Cbox,SHub] =mrQ_boxScaleGlobLinear(ScaleMat)
%[Cbox,SHub] =mrQ_boxScaleGlobLinear(ScaleMat)
%we will find the scale of each box so it will agreee with the values in the
%overlap boxes ( a scalar term). We assume that the boxPD  are free from gain and the only different is a scaler that need to be fit.
% all the diffent scalr or each two box are input ScaleMat.
% boxII= boxJJ*Ratio
 %               scaleFactor(ii,jj)=1./Ratio;
  %              scaleFactor(jj,ii)=Ratio;
                
%we will right a minimization eqation so all the scalares thogther will be adjasted.
%the idea is  box(i) -box(j)*scaler equal 0. we will build such set of eqations for all
%the overlap boxes and solve them thogther in the end.
% one box will be a refernce box.
%in the case of the reference we deside that scaler is 1;
%
% simulation of the problem we solve is at SimLinBoxJoin.m 
%
%
% AM (C) Stanford University, VISTA

%% find the bigest network of connected boxes

% what is the tipical ratio between two box?
MScale=median(ScaleMat(ScaleMat>0));
% if a ratio is more then 2 fold the usual it might be wrong. it might be
% better to get ride of those box. tipically this is only a small fraction
% then the box >>1%. this operation is usful for the solving the linear
% system without outlayers.
ScaleMat(ScaleMat>MScale*2)=0;
ScaleMat(ScaleMat<MScale*0.5)=0;



 [S, C]= graphconncomp(double(sparse(logical(ScaleMat))));
 for ii=1:S
 NudeN(ii)= length(find(C==ii));
 end
  [N,ind]=sort(NudeN);
  
% all the box that are conectet in the bigest network
boxT0Use=(C==ind(end)) ;

% build the  matrix for the scale linear calculation
LinScaleMat=zeros(size(ScaleMat));

% calcualte how many conction each box got
NConect=sum(logical(ScaleMat),1);
 Hub= (NConect==max(NConect)) & boxT0Use;

 %find the nude that highly conected and select one
 SHub=find(Hub);SHub=SHub(1,1);
LinScaleMat(SHub,SHub)=1;

%% build the eqations matrix
% say refernce=1; and there are 3 boxs
%
% MAT=
%
%             box(reference)                                      =1
%          2* box(reference) -scale(2,1)*box(2) -scale(3,1)*box(3)=0
% -scale(1,2)*box(reference)          +2*box(2) -scale(3,2)*box(3)=0
% -scale(1,3)*box(reference) -scale(2,3)*box(2)          +2*box(3)=0


for   ii=find(boxT0Use)
    if (ii ~=SHub)
       
 %           boxII= boxJJ*Ratio
 %               scaleFactor(ii,jj)=1./Ratio;
 %              scaleFactor(jj,ii)=Ratio;

     
 LinScaleMat(ii,:)=-ScaleMat(ii,:);
 LinScaleMat(ii,ii)=NConect(ii);
    end
end


BoxLocation=find(boxT0Use);

y=zeros(size(LinScaleMat,1),1);
y(SHub)=1;
y=y(BoxLocation);
Mat=LinScaleMat(BoxLocation,BoxLocation);

%solve it as multi linear eqation
C=pinv(Mat'*Mat)*Mat'*y;
% the C we need is one over the one we fit


Cbox=zeros(size(LinScaleMat,1),1);

Cbox(BoxLocation)=C;



