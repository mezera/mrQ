function [nopossible]=errlocalGainUC_vQualety(x,box,Poly,coefdat,use,coils,boxS)
%
%[nopossible]=errlocalGainUC_vQualety(x,box,Poly,coefdat,use,coils,boxS)
% Input:
% x     -the fitted parameters - the coil polynomials coefitiones
%
% box     - the raw images of the different coils before the fit
%
% Poly   - the polynomial we used to fit (the one we multipied with x )
%
%coefdat - the different coils images coraltion within the box
%
%coils   - the number of coils
% 
%boxS   - the imaging box dimention where the data sit


% # we will check if the result is relevant and try to find the problematic
%coils
%we will run 4 checks
% 1.we will check if the coil max is in the edge of the box.
% 2. we will check if the fit of one of the box is totly off (more then
% double the median std of the error
% 3. check that there are no nan or inf in the result (can be a result of
% bad input coil or bad Gain fit
% 4. check that the Gain is not more coralate then the data.
%if the coil found max in

%a list of the coils that have a fit  that is not possible
nopossible=zeros(coils,1);

%%
%check 1: is the max in the edge of the box. if it is not it must be wrong
%(see below)
for i=1:coils
 Gain1(:,:,:,i) = reshape(Poly*x(i,:)',boxS);
end
%nested function (see below)
nopossible=check_possible(Gain1);

clear Gain1

%%
%check 2 the error std is not an outlayer
% The gain for each of the coils, estimated by the parameter x

%calculate the gain of each coil
Gain = Poly*x(:,:)';
%mreshape it from vector to 3D images (4D) 3Dxcoils
Gain =reshape(Gain(use),[],coils);
%calcualte the PD of each box 
Val = box./Gain;% 4D pd of each coil

% only the one that pass the chaeack above
Val=Val(:,nopossible==0);
coils=length(find(nopossible==0));
%the PD of all should be the same 
%we mesure the error
err1=(((Val - repmat(mean(Val,2),1,coils)) )./repmat(mean(Val,2),1,coils));

%cheack the mean and std
stds=std(err1);
Mstd=median(stds);
% is a coil is off then this is a bad fit
if find(stds>Mstd.*2);nopossible(find(stds>Mstd.*2))=1;end
%%
% check 3 is there  nan or inf
if  find(isnan(Val))  
[d y1]=ind2sub(size(Val),find(isnan(Val)));
nopossible(y)=1;
end
if  find(isinf(Val))  
[d y1]=ind2sub(size(Val),find(isinf(Val)));
nopossible(y)=1;
end

%%
%check 4 gain coralation between coil and data
% we fit the gain under the condtion that after correction the coil bias is
% always less coralate then when each of them multiple by the same PD.
%if this is not the case then the fit is bad and some bias was added by our
%fit

%calculate the fitted coil gaincoraltion
coefG =  tril(corrcoef(Gain),-1);

% Outside of this routine we already calculated
%   coefdat =tril(corrcoef(box),-1);
% Here we compare the two correlation coefficients
% The correlation coefficients of the data (coefdat) should be larger than
% the corr coef of the gain (coefG).

err2 = (coefdat-coefG)./abs(coefdat);
% we will go over the coil 
%If the difference is positive, there is no problem.
%if not we need to find the problematic coil and deside which one to excude 

err2(err2>0)=0;

if find(err2<-0.05);
   [C1(1,:) C2(1,:)]= ind2sub(size(err2),find(err2<-0.01));
   nogood=find(nopossible);
       if length(nogood)~=0
   for i=1:length(nogood)
     
       C2(find(C2==nogood(i)))=0;
   C1(find(C1==nogood(i)))=0;
   end
       end
   
   while find(C1.*C2)
       tt1=[C1(find(C1.*C2)) C2(find(C1.*C2))];  
    tt=sort(unique(tt1),'descend');  tt=tt(tt>0) ;
    if length(tt)>0
        clear count
    for i=1:length(tt)
        count(i)=length(find(tt1==tt(i)));
    end
    [d I]=max(count);
    end
   nopossible(tt(I))=1; 
    C1(find(C1==tt(I)))=0;
        C2(find(C2==tt(I)))=0;
end
end

%%
%nested function 
%all the boxes are away from the coil.
%therefor a fitted  coil gain that is minimum or maxsimum value is not on
%the edge of the fitted box must be wrong. it is like the fitted box was
%pysicaly inside the coil ....
function [notposible]=check_possible(Gain)
sz=size(Gain);
notposible=zeros(sz(4),1);

for i=1:sz(4),
    tmp=Gain(:,:,:,i);
    [cx cy cz]=ind2sub([sz(1) sz(2) sz(3)],find(tmp==max(tmp(:))));

    if (min(cx)==1 || max(cx)==sz(1) || min(cy)==1 || max(cy)==sz(2)  || min(cz)==1 || max(cz)==sz(3))
        notposible(i)=0;
    else
        notposible(i)=1;
    end;
end;



%% for debuging and visualization of the fit and the data from vector to 3D and 4D

% %%
%  ffVal=zeros([boxS coils]); ffbox=zeros([boxS coils]);
%  bm=zeros(boxS);bm(use(:,1))=1;bm=logical(bm);
%  for i=1:coils
%      tmp=zeros(boxS);
%      tmp(bm)=Val(:,i);
%      ffVal(:,:,:,i)=tmp;
%     tmp=zeros(boxS);
%      tmp(bm)=box(:,i);
%      ffbox(:,:,:,i)=tmp; 
% end
%           tmp(bm)=mean(Val,2);
%     
