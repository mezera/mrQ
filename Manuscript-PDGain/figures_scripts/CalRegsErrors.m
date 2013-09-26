%% 
% This script call for the simulate and fitting presider many time fit
% different random noise  PD distributions, and coils functions
% calculate the error of each of them.
tic
for ii=1:25
    for loc=1:5
        for PD=1:10
      [ Err_cor(ii,loc,PD,:),   Err_ridge(ii,loc,PD,:), Err_T1reg(ii,loc,PD,:) , Err_Noreg(ii,loc,PD,:)]=SimAllmodelFit(2,loc,num2str(PD-1)) ;     
        end
    end
end
toc;
save('/home/avivm/Documents/PD_article/figures/fig3MultiRegFits','Err_cor','Err_ridge','Err_T1reg','Err_Noreg')