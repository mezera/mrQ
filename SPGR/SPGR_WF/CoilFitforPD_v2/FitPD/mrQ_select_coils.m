function coils=mrQ_select_coils(nUseCoils,MaxcoilNum,M0_v)
% check the correlation btween the coils and select the nUseCoils form the
% 1:MaxcoilNum in M0_v mo matrix (voxel X coils). 

% We use this algorithm to
c = nchoosek(1:MaxcoilNum,nUseCoils);          % All the potential combinations of nCoils
Cor = ones(nUseCoils,size(c,1))*100;   % Initiate the Cor to max
for kk=1:size(c,1)             % loop over coils combinations
    
    % Correlations between the the measured M0 for each of the coils in
    % this combination.
    A = (corrcoef(M0_v(:,c(kk,:)))); 
    
    % Sum the abs correlation after correcting for number of coils and the
    % identity correlations
    Cor(nUseCoils,kk) = sum(sum(abs(triu(A) - eye(nUseCoils,nUseCoils))));
    
end

% Find the minimum correlation (min abs corr give us the set with corr that
% are closer to zero).  Choose those coils.
[v, minCor] = min(Cor(nUseCoils,:)); 
coils       = c(minCor,:);


