function PDinit=Get_PDinit(T1reg,R1,Guesstype)
%PDinit=Get_PDinit(lambda1,T1,Guesstype)
% estimate a PD from the data
%  Guesstype=1 or T1reg=0 use the sum of sqr (or M0) as PD
%  Guesstype=2 randum start
% Guesstype=3 start or T1reg>0  with the PD T1 relation as describe in the literature.
%
% AM  & BW VISTASOFT Team, 2013

if notDefined('Guesstype')
Guesstype=0;
end

if (T1reg==0 || Guesstype==1) % sum of squr
PDinit = sqrt(sum(M0.^2,2));  
end

if  Guesstype==2 %rand PD
PDinit = rand(size(M0,1),1);
end

if (T1reg>0 || Guesstype==3)% T1 PD linear relations
    
          PDinit=1./(R1*0.42+0.95); %this is the teortical T1 PD relationship see reference at Mezer et. al 2013
end

