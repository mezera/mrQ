%% 
addpath(genpath(fullfile(mrqRootPath)));

%%
% This script call for the simulate and fitting presider many time fit
% different random noise  PD distributions, and coils functions
% calculate the error of each of them.

%
tic
for ii=1:25
    for loc=1:5
        for PD=1:10
      [ Err_cor(ii,loc,PD,:),   Err_ridge(ii,loc,PD,:), Err_T1reg(ii,loc,PD,:) , Err_Noreg(ii,loc,PD,:)]=SimAllmodelFit(2,loc,num2str(PD-1)) ;     
        end
    end
end
toc;
save('/home/avivm/Documents/PD_article/figures/fig3/MultiRegFits1','Err_cor','Err_ridge','Err_T1reg','Err_Noreg')


%
% Elapsed time is 33781.445993 seconds.
%
% avrage along coils location and random noise. to each PD case
Cor_2=Err_cor(:,:,:,2);C(:,2)=squeeze(median(median(Cor_2)))';
Rig_2=Err_ridge(:,:,:,2);R(:,2)=squeeze(median(median(Rig_2)))';
T1_2=Err_T1reg(:,:,:,2);T(:,2)=squeeze(median(median(T1_2)))';
No_2=Err_Noreg(:,:,:,2);N(:,2)=squeeze(median(median(No_2)))';

Cor_1=Err_cor(:,:,:,1);C(:,1)=squeeze(median(median(Cor_1)))';
Rig_1=Err_ridge(:,:,:,1);R(:,1)=squeeze(median(median(Rig_1)))';
T1_1=Err_T1reg(:,:,:,1);T(:,1)=squeeze(median(median(T1_1)))';
No_1=Err_Noreg(:,:,:,1);N(:,1)=squeeze(median(median(No_1)))';
