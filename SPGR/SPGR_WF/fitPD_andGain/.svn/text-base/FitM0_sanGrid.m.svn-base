function FitM0_sanGrid(opt,jumpindex,jobindex);

j=0;
st=1 +(jobindex-1)*jumpindex;
ed=st+jumpindex-1;
if ed>length(opt{1}.wh), ed=length(opt{1}.wh);end;
options = optimset('Algorithm', 'levenberg-marquardt','MaxIter',100,'MaxFunEvals',inf,'Tolx',1e-12);

%options =  optimset('LevenbergMarquardt','on','Tolx',1e-12,'TolF',1e-12,'MaxIter',100,'MaxFunEvals',inf);
nn=1;
for i= st:ed,
it=1;exitflag=0;x0=opt{i}.x0;
resnorm(nn)=0;resnorm0=0;
        while (exitflag<1 && it<10 && resnorm0/resnorm(nn)~=1 )
resnorm0=resnorm(nn);
[res(:,:,nn), resnorm(nn),dd,exitflag] = lsqnonlin(@(par) errlocalGainUC(par,opt{i}.inDat,opt{1}.boxS,opt{1}.Poly,opt{1}.inNan,opt{1}.numIn),x0,opt{1}.lb,opt{1}.ub,options);
%[res(:,:,nn), resnorm(nn),dd,exitflag] = lsqnonlin(@(par) errlocalGainCSF(par,opt{i}.inDat,opt{1}.boxS,opt{1}.Poly,opt{1}.inNan,opt{1}.numIn,opt{i}.known),x0,opt{1}.lb,opt{1}.ub,options);
it=it+1;
x0=res(:,:,nn);
  end
nn=nn+1;


end;
name=[opt{1}.outDir opt{1}.name '_' num2str(st) '_' num2str(ed)];
save(name,'res','resnorm','st','ed')



    