function [csfnm]= mrQM0fit_isBoxcsf(opt,brainMask,fb,csf) 


sz=size(brainMask);


[Xx Yy Zz]=MrQPD_boxloc(opt,fb);


if (sz(1)< Xx(1) || sz(1)< Xx(2)) 
    box=[]; 
elseif (sz(2)< Yy(1) || sz(2)< Yy(2))    
box=[]; 
elseif (sz(3)< Zz(1) || sz(3)< Zz(2))    
box=[]; 
else
    
box=csf(Xx(1):Xx(2),Yy(1):Yy(2),Zz(1):Zz(2));

end;
csfnm=length(find(box));
