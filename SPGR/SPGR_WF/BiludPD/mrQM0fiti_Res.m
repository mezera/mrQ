function [M0f,donemask,Avmap,M0ft,M0ftfull] =mrQM0fiti_Res(opt,fb,res,inDat,M0f,donemask,known,Avmap,M0ft,M0ftfull,resnorm,exitflags)

 rel=find(known(2,:)<120 & known(2,:)>70 );
    
    if length(find(rel))<10;
    return
    end
[Xx Yy Zz]=MrQPD_boxloc(opt,fb);

for i=1:size(inDat,4),%opt.numIn
    %
    Gain1(:,:,:,i) = reshape(opt.Poly*res(i,:)',opt.boxS);
    Val(:,:,:,i) = inDat(:,:,:,i)./Gain1(:,:,:,i);
end;
%[inDat,Val,skip]=mrQcheck_possible(opt,Gain1,fb,inDat,Val);

% if skip==1
%     return;
% end;

if(exist('known','var'))
    ResVal=median(Val,4); 
    in=ResVal(known(1,rel))./known(2,rel);
    med=median(in);rel=find(in<med*5 & in>med/5);
    [d dd]= ksdensity(in(rel), linspace(min(in(rel)),max(in(rel)),5000));%[min(in(rel)):0.01:max(in(rel))] );
    dc=dd(find(d==max(d)));

    
    % dc=mean(ResVal(known(1,:))./known(2,:));
    
    ResVal=median(Val,4)./dc;
else
    ResVal=median(Val,4);
end
    box=M0ft(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
     
    if(exist('Avmap','var'))
        Av=Avmap(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
        

    else
        Av=ones(size(box));
    end
    
    
%     if mean(ResVal(find(ResVal)))>105
%         keyboard;
%     end;
    
    %box1=box;
    if (find(ResVal<0  | ResVal>170))

                      %  keyboard;

                      ermin=find(ResVal<0);
                      ermax=find(ResVal>200);
                      %tt=find(box(ermax)==0);
                      %box(ermax(tt))=200;
                      
                      ResVal(ermax)=box(ermax);
                      
                      ResVal(ermin)=box(ermin);
                      
        Av(ermin)=Av(ermin)-1;
        Av(ermax)=Av(ermax)-1;


        if find(Av<0);
            Av(find(Av<0))=0;
        end
        
    end;
    box=(box.*Av+ResVal)./(Av+1);
    
    
    % M0ftfull(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3),ssz1)=ResVal;
   
    
    % showMontage(box); caxis([0 1.2])
    % showMontage(box1); caxis([0 1.2])
    
   % M0ft(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3))=box;
    M0ft(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=box;
    %M0f(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3))=box;
    M0f(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=box;
    
    Avmap(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=Avmap(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))+1;

%     Avmap(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3))=...
%         Avmap(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3))+1;
%     
%    ssz1=max(Avmap(:)); 

    
    donemask(fb(1),fb(2),fb(3))=1;
   if find(isnan(M0f))
    haa=1;
   end

  
   
    ssz=size(M0f);
 ssz1=max(Avmap(:)); %max(max(Avmap(:)),ssz1);
  ssz(4)=size(M0ftfull,4);
%  loc=Avmap(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3));
  loc=Avmap(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
% [Xx,Yy,Zz] = ndgrid(X(fb(1),fb(2),fb(3))-HboxS(1):X(fb(1),fb(2),fb(3))+HboxS(1),...
 %                      Y(fb(1),fb(2),fb(3))-HboxS(2):Y(fb(1),fb(2),fb(3))+HboxS(2),...
  %                     Z(fb(1),fb(2),fb(3))-HboxS(3):Z(fb(1),fb(2),fb(3))+HboxS(3));
 
                   [Xx1,Yy1,Zz1] = ndgrid(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
 
     wh=  sub2ind(ssz,Xx1(:),Yy1(:),Zz1(:),loc(:));
  
  
 if ssz(4)<ssz1
 tmp=zeros(ssz);
 tmp(:,:,:,1:ssz1)=M0ftfull; 
 tmp(wh)=ResVal(:)';  
 M0ftfull=tmp;
 
 else
  M0ftfull(wh)=ResVal(:);    
 end;
 %



     
     
     
     

               