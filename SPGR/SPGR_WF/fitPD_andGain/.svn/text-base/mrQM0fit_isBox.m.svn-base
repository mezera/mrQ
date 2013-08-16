function [empty,Avmap] = mrQM0fit_isBox(opt,brainMask,fb,Avmap)
% 
% [empty,Avmap] = mrQM0fit_isBox(opt,brainMask,fb,Avmap)
% 
% Comments to come...
% 
% INPUTS:
%       opt     - 
%    brainMask  - 
%       fb      - 
%     Avmap     - 
% 
% OUTPUTS:
%       empty   - 
%       Avmap   -
% 
%      
% 
% (C) Stanford University, VISTA
% 

%%

sz = size(brainMask);

[Xx Yy Zz] = MrQPD_boxloc(opt,fb);


if (sz(1)< Xx(1) || sz(1)< Xx(2))
    box=[];
elseif (sz(2)< Yy(1) || sz(2)< Yy(2))
    box=[];
elseif (sz(3)< Zz(1) || sz(3)< Zz(2))
    box=[];
else
    box=brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
end

% if (60> Xx(1) &&  60< Xx(2)  && 176> Yy(1) &&  176< Yy(2) && 100> Zz(1) &&  100< Zz(2)   )
%     c=1;
% end

cutT=length(box(:)).*0.05; % we will work only with boxes that have data on at least 5% of it. and at least 200 voxels

if (length(find(box))>cutT && length(find(box))>200)
    empty=0;
    Avmap(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))=Avmap(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2))+1;
    
elseif length(find(box))<1
    empty=1;
else
    empty=-1;
    
end