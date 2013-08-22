load('X_validationFit_1.mat')


for Loc=1:sampleLocation
    
    for NS=1:length(nSamples)
        
        for C=LastCoil
                                        MaxkFold=C;

            for K=2:MaxkFold
PDRMSE=RMSEAll(:,Loc,NS,C,K);
                X_valdationErrF= X_valdationErrAll(:,:,Loc,NS,C,K);
                mrvNewGraphWin;
                subplot(1,2,1)

                 semilogy(lambda1,X_valdationErrF(1,:),'*-')
                best = find(X_valdationErrF(1,:)==min(X_valdationErrF(1,:)));
                  ylabel('sum of abs err')
                xlabel('lambda')
              title(  ['Choosing lambda is ' num2str(lambda1(best))  ' gave PDerr ' num2str(PDRMSE(1) ) ' for  Ncoils: '  num2str(C)  'and Kfold '  num2str(K) 'for  Nsample:' num2str(nSamples(NS))   ]) 
                subplot(1,2,2)
         semilogy(lambda1,X_valdationErrF(2,:),'*-')
                best = find(X_valdationErrF(2,:)==min(X_valdationErrF(2,:)));
                              title(  ['Choosing lambda is  ' num2str(lambda1(best))  ' gave PDerr ' num2str(PDRMSE(2) )]) 

                ylabel('RMSE')
                xlabel('lambda')
keyboard
close
            end
        end
    end
end




mrvNewGraphWin;
subplot(1,2,1)
imagesc(( squeeze(RunTime(1,1,:,:))))
ylabel('coils')
xlabel('fold')
title('time in sec for 3 samples')
colorbar
subplot(1,2,2)
imagesc(( squeeze(RunTime(1,:,:,2))))
ylabel('Nsample')
xlabel('coils')
title('time in sec for 2flod  X validation')
colorbar


mrvNewGraphWin;
subplot(1,2,1)
imagesc(( squeeze(RMSEAll(2,1,3,:,:))))
ylabel('coils')
xlabel('fold')
title('PD err for 4 samples')
colorbar
 caxis([0 4])
subplot(1,2,2)
imagesc(( squeeze(RMSEAll(2,1,:,:,2))))
ylabel('Nsample')
xlabel('coils')
title('PD err for2 flod  X validation')
colorbar
 caxis([0 4])
