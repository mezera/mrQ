

MM=M;
MM(isnan(MM))=0;

coefdat=tril(corrcoef(MM),-1);

MM=MM>0;

WHM=find(mask);


for i=1:coils
    [Py{i}.Poly,str] = constructpolynomialmatrix3d(size(BM),WHM(MM(:,i)),degrees);
end



for i=1:coils
   in=zeros(size(BM));
    in(mask)=M(:,i);
    inMask=zeros(size(BM));
    inMask(WHM)=1;
    inMask=inMask & in>0;
    [params1,gains,rs] = fit3dpolynomialmodel(in,inMask,degrees);
    X0(:,i)=params1';
    
    
end
    
   Gain=nan(length(find(mask)),coils);
   
 a=version('-date');
if str2num(a(end-3:end))==2012
    options = optimset('Algorithm', 'levenberg-marquardt','Display', 'iter' ,'Tolx',1e-6,'MaxIter',100,'MaxFunEvals',inf);
else
    options =  optimset('LevenbergMarquardt','on','Display', 'iter','Tolx',1e-6,'MaxIter',100,'MaxFunEvals',inf);%'TolF',1e-12
    
end

c
%err=errAllcoils(x,Gain,a,M,MM,coefdat)
[res, resnorm,dd1,exitflag] = lsqnonlin(@(par) errAllcoils(par,Gain,Py,M,MM,coefdat,coils),double(X0),[],[],options);

