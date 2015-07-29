function [empty] = mrQ_isDataBox(opt,brainMask,fb,Inclusion_Criteria)
% [empty] = mrQ_isDataBox(opt,brainMask,fb,Inclusion_Criteria)
% 
% Comments to come...
% 
% ~INPUTS~
%                     opt:  Structure of parameters.
%               brainMask: 
%                      fb: 
%      Inclusion_Criteria:  Default is [0.05 200].
% 
% ~OUTPUTS~
%                   empty:  
%      
% 
% (C) AM Stanford University, VISTA
% 

%%

sz = size(brainMask);

Xx(1)=opt.X(fb(1),fb(2),fb(3))-opt.HboxS(1);
Xx(2)=opt.X(fb(1),fb(2),fb(3))+opt.HboxS(1);
Yy(1)=opt.Y(fb(1),fb(2),fb(3))-opt.HboxS(2);
Yy(2)=opt.Y(fb(1),fb(2),fb(3))+opt.HboxS(2);
Zz(1)=opt.Z(fb(1),fb(2),fb(3))-opt.HboxS(3);
Zz(2)=opt.Z(fb(1),fb(2),fb(3))+opt.HboxS(3);

%check this is not outside the image
if (sz(1)< Xx(1) || sz(1)< Xx(2))
    box=[];
elseif (sz(2)< Yy(1) || sz(2)< Yy(2))
    box=[];
elseif (sz(3)< Zz(1) || sz(3)< Zz(2))
    box=[];
else
    box=brainMask(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));
end

if notDefined('Inclusion_Criteria')
    Inclusion_Criteria=[0.05 200];
end

% if (60> Xx(1) &&  60< Xx(2)  && 176> Yy(1) &&  176< Yy(2) && 100> Zz(1) &&  100< Zz(2)   )
%     c=1;
% end

% We will work only with boxes that have data on at least 5% of it, and at least 200 voxels
cutT=length(box(:)).*0.05; 

% Check it's not almost empty from mask voxel
if (length(find(box))>=length(box(:)).*Inclusion_Criteria(1) && length(find(box))>Inclusion_Criteria(2))
    empty=0;
elseif length(find(box))<1
    empty=1;
else
    empty=-1;
end