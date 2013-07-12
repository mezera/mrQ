function PD= pdGuessFromM0(M0,GussTypeFlug)

if notDefined('GussTypeFlug')
GussTypeFlug='sumofsqr')
end

nCoils=size(M0_v,2);
nvoxel=size(M0_v,1);

 % This is the PD that will be returned if we won't have multi coil
    % information
    if ischar(GussTypeFlug)
        GussTypeFlug = mrvParamFormat(GussTypeFlug);
    end
        switch GussTypeFlug
            case {'sumofsqr', '1'}
   
    PD1 = sqrt(sum(M0_v.^2,2));
    PD1= PD1./ mean(PD1);   % There are options - e.g., we could set PD(1) to 1
    
     case {'demeansum', '2'}
    Mn=mean(M0_v);
         M0_=M0_v./repmat(Mn,nvoxel,1);
         PD2=mean(M0_,2);
  
         
     case {'demeansumUC', '3'}
         
         
         
        end