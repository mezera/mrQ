function PD= pdGuessFromM0(M0,GussTypeFlug)

if notDefined('GussTypeFlug')
GussTypeFlug='sumofsqr')
end

nCoils=size(M0,2);
nvoxel=size(M0,1);

 % This is the PD that will be returned if we won't have multi coil
    % information
    if ischar(GussTypeFlug)
        GussTypeFlug = mrvParamFormat(GussTypeFlug);
    end
        switch GussTypeFlug
            case {'sumofsqr', '1'}
   
    PD = sqrt(sum(M0_v.^2,2));
    PD = PD./ mean(PD);   % There are options - e.g., we could set PD(1) to 1
    
     case {'demeansum', '2'}
    Mn=mean(M0);
         M0=M0./repmat(M,nvoxel,1);
         PD=mean(M0,2);
  
         
     case {'demeansumUC', '3'}
         
         
         
        end