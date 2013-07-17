function [S S1]=RatioVar(M)
SZ=size(M);

for i=1:SZ(2)
    for j=1:SZ(2)

S(i,j)=(std( M(:,i)./M(:,j)) +  std( M(:,j)./M(:,i))) /2;
S1(i,j)=std(M(:,i)) - std(M(:,j))   ;

    end
end

