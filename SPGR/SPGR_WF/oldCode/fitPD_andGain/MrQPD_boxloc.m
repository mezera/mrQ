function [Xx Yy Zz,skip] = MrQPD_boxloc(opt,fb)
% 
% [Xx Yy Zz,skip] = MrQPD_boxloc(opt,fb)
% 
% Comments to come...
% 
% INPUTS:
%       opt - 
%       fb  - 
% 
% OUTPUTS:
%       Xx  - 
%       Yy  -
%       Zz  - 
%     Skip  - 
% 
%      
% 
% (C) Stanford University, VISTA
% 
% 

%%

skip=0;

if opt.donemask(fb(1),fb(2),fb(3))==-1000
    Xx=0;Yy=0; Zz=0;
    
    display ('this is an empty box (no brain mask)')
    skip=1; return
end

Xx(1)=opt.X(fb(1),fb(2),fb(3))-opt.HboxS(1);
Xx(2)=opt.X(fb(1),fb(2),fb(3))+opt.HboxS(1);
Yy(1)=opt.Y(fb(1),fb(2),fb(3))-opt.HboxS(2);
Yy(2)=opt.Y(fb(1),fb(2),fb(3))+opt.HboxS(2);
Zz(1)=opt.Z(fb(1),fb(2),fb(3))-opt.HboxS(3);
Zz(2)=opt.Z(fb(1),fb(2),fb(3))+opt.HboxS(3);


if opt.donemask(fb(1),fb(2),fb(3))==-2000
    out=0;
    yes=zeros(6,1);
    sz=size(opt.donemask);
    jump(1)=Xx(2)-Xx(1);
    jump(2)=Yy(2)-Yy(1);
    jump(3)=Zz(2)-Zz(1);
    jump=round(jump./1.2);
    if  (fb(1)>1 && out==0);
        if opt.donemask(fb(1)-1,fb(2),fb(3))==0;
            Yy(1)=Yy(1)-jump(1); Yy(2)=Yy(2)-jump(1);
            out=1;
        end;
    end;
    if  (fb(1)<sz(1) && out==0);
        if opt.donemask(fb(1)+1,fb(2),fb(3))==0;
            Yy(1)=Yy(1)+jump(1); Yy(2)=Yy(2)+jump(1);
            out=1;
        end;
    end;
    if  (fb(2)>1 && out==0);
        if opt.donemask(fb(1),fb(2)-1,fb(3))==0;
            Xx(1)=Xx(1)-jump(2); Xx(2)=Xx(2)-jump(2);
            out=1;
        end;
    end;
    if  (fb(1)<sz(2) && out==0);
        if opt.donemask(fb(1),fb(2)+1,fb(3))==0;
            Xx(1)=Xx(1)+jump(2); Xx(2)=Xx(2)+jump(2);
            out=1;
        end;
    end;
    if  (fb(3)>1 && out==0);
        if opt.donemask(fb(1),fb(2),fb(3)-1)==0;
            Zz(1)=Zz(1)-jump(3); Zz(2)=Zz(2)-jump(3);
            out=1;
        end;
    end;
    if  (fb(3)<sz(3) && out==0);
        if opt.donemask(fb(1),fb(2),fb(3)+1)==0;
            Zz(1)=Zz(1)+jump(3); Zz(2)=Zz(2)+jump(3);
            out=1;
        end;
    end;
if out==0, % This is a box that is an island lets just sikp it
skip=1;
end

end;
