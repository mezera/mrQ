function OutPut = pdBiLinearFit_2(M0_v,pBasis,Lambda,maxLoops,sCriterion,PD,plotFlag,TruePar,D,HoldforCV,W)
% OutPut = pdBiLinearFit(M0_v,pBasis,Lambda,maxLoops,sCriterion,PD,plotFlag,TruePar)
%
% This fancrtion running a while loop solving the M0=GPD bi-linear problem.
% using ridge regression with regularisation coeffisent Lamdea for the gain
% problem. im each interval of the loop. the coil gain are estimated given
% a PD values. then the PD is evaluated again given the coil gain sulotion.
% the loop stop when the solotion convarge up to the stopping criterion or
% when it reach the maxsimal iteration
%
% Input
%  M0_v:         Columns containing the 3D M0 values in a vector (nVoxels x nCoils)
%  pBasis:       Polynomial basis (nPositions x nCoef)
%  Lambda        A ridge regration regularization term
%  maxLoops     maxsimal iteration of the while loop
%  sCriterion    Stopping criterion to stop the while loop (this  define what when the problem converge)
%  PD                  initial PD sulotion
%  plotFlag         When true (1) output a set of images showing the progess of the fit
%  TruePar       the true coil gain coepicent parameters if they are
%                       known. like in the casses of   simulation or phantoms.
% D                    A wighted identity matrix that will wighting  the ridge
%                       regration for each parameters D=(NcoefXnCoef)   (Ncoef=size(pBasis,2)); .
%                       D is a diagonal matrix position 1,1, is wight on parameter 1 position 2,2
%                       in a wight on parameter 2 ec. position (n,n) wight on parameter n.
%HoldforCV     The fraction of the data that will be hold fro cross.HoldforCV isa number between 0 to 1 HoldforCV>0 HoldforCV<1
%W                     A diagonal matrix  (nVoxels X nVoxels)  of the wight
%                         of each raw im M0 for the lsq fit
%
% OutPut is a stracture with a list of outputs:
%  OutPut.PD             -  the final PD values
%  OutPut.Gn             -  the final coils gain values
%  OutPut.g                -  the final coils gain coeficents  (G = pBasis*g)
%  OutPut.PDchange -  a vector of the PD change in each step
%  OutPut.NumOfIter - number of iteration
%  OutPut.convarge    - 1 if convarge 0 if not
%
%   AM/BW Copyright VISTASOFT Team 2013

%% intiate parameters

if notDefined('maxLoops'), maxLoops = 100; end
if notDefined('sCriterion')
    sCriterion = 1e-4;  % Stopping criterion
end
if notDefined('HoldforCV') HoldforCV=0;end


nCoils=size(M0_v,2);
nvoxel=size(M0_v,1);
nPolyCoef = size(pBasis,2);




if HoldforCV>0  && HoldforCV<1
    % cross validation we selecrt a random set of mesurament fot the fit
    % and a random set to fit on
     voxelLocations=randperm(nvoxel);
     cut=round(nvoxel*HoldforCV);
     
     CVvox=voxelLocations(1:cut);
     Fitvox=voxelLocations(cut+1:end);
     nvoxelf=length(Fitvox);

else
    % no cross validation we will use all the mesrament to fit and to
    % mesure the godness of the fit
        nvoxelf=nvoxel;
    Fitvox=1:nvoxel;
    CVvox=Fitvox;
    end
    

if notDefined('plotFlag'), plotFlag = 0; end


if notDefined('D'),  D=eye(nPolyCoef); end
if notDefined('W'),  W=[]; end



% loop and solve by ridge regration
k = 0;                % number of iteration
tryagain = 1;         % go for the while loop
PDchange = zeros(1,maxLoops);
M0change = zeros(1,maxLoops);
M0Fit = zeros(1,maxLoops);
M0CVFit = zeros(1,maxLoops);  
stablePD=0;
stableM0=0;

convarge=0;
%% Initalize solution
if notDefined('PD')
    % This is the PD that will be returned if we won't have multi coil
    % information
   
g0 = zeros(nPolyCoef,nCoils);
for ii=1:nCoils
    g0(:,ii) = pBasis(Fitvox,:) \ M0_v(Fitvox,ii);  % Polynomial approximation
end
[PD, G] = pdEstimate(M0_v, pBasis, g0);
end

% You must first divide G by PD(1).  If you divided PD first, then you
% change PD and you no longer are using the valid PD(1).
G  = G  .* PD(1);
PD = PD ./ PD(1);

M0 = G.*repmat( PD(:),1,nCoils);
%% This plots the initial condition
if plotFlag==1
    figH = mrvNewGraphWin;
    if notDefined('TruePar')
        Par=g0;
        GetPar=1;
        str='last interval coefisents';
    else
        Par=TruePar;
        str='True gains';
        GetPar=0;
    end
    CoefNorm=Par(1)./g0(1);
    
    % Plot the starting point
    figure(figH);
    set(figH,'Name', ['Loop ' num2str(k)]);
    subplot(1,2,1); plot(g0(:).*CoefNorm,Par(:),'.')
    identityLine(gca);
    xlabel('Estimated gains'); ylabel(str);
    
    subplot(1,2,2); plot(PD,'.');set(gca,'ylim',[min(PD(:))*0.9 max(PD(:))*1.1]);
    ylabel('Estimated PD');xlabel('voxel');
    
end

%% This is the bilinear alternating solution
while tryagain == 1
    k = k+1; % count Ridge regression steps
    
    % fit linear equation coefficients using ridge regression
    g = RidgeRegressCoilfit(PD(Fitvox), Lambda, M0_v(Fitvox,:), pBasis(Fitvox,:),D,W);
    
    % Calculate the coil gain of each coil over space from the estimate
    % Calculate PD for each coil
%    [PDn, Gn] = pdEstimate(M0_v(Fitvox,:), pBasis(Fitvox,:), g);
        [PDn, Gn] = pdEstimate(M0_v(:,:), pBasis(:,:), g);  %as we solve voxel by voxel we can run on all of them

    
    Gn  = Gn  .* PDn(1);
    PDn = PDn ./ PDn(1);
    
    M0n = Gn.*repmat( PDn,1,nCoils);
    % Check if the new estimate differs from the one before or it's
    % converged
    PDchange(k) = std(PD - PDn);
    M0change(k) = std(M0(:) - M0n(:));
    
    CVterm=M0_v(CVvox,:)-M0n(CVvox,:);
      M0CVFit  (k) = std(CVterm(:));
 FitTerm=M0_v(Fitvox,:)-M0n(Fitvox,:);
        M0Fit(k) = std(FitTerm(:));

    %if PDchange(k) < sCriterion;
    if M0Fit(k) < sCriterion   || min(M0CVFit(1:k)) -M0CVFit(k) < -sCriterion 
        % if the data is fitted stop
        
        % We could check the gains, rather than PD, or both
        % if std(G-Gn)<0.01
        
        % if the two solutions are stable --> stop
        tryagain=0;
        convarge=1;
                            fprintf(' rich the minimum for M0 fit CV  ');

    else
        % Keep going.
        
        % Update the new PD and and estimated M0
        PD = PDn;
        M0 = M0n;
        
        if plotFlag==1
            %keyboard
            CoefNorm=Par(1)./g(1);
            % plot the new estimations for coil gain and PD
            figure(figH);
            set(figH,'Name',[ 'Loop ' num2str(k) ] );
            subplot(1,2,1);plot(g(:).*CoefNorm,Par(:),'.')
            identityLine(gca);
            xlabel('Estimated gains'); ylabel(str);
            
            subplot(1,2,2);
            plot(PDn,'.'); set(gca,'ylim',[min(PDn(:))*0.9 max(PDn(:))*1.1]);
            ylabel('Estimated PD');xlabel('voxel');
            if  GetPar==1
                Par = g;
            end
        end
    end
    
%      if k>10 &&  (max(PDchange(k-10:k))-PDchange(k))./PDchange(k)<0.01
%             % If PD stable to within 1 percent, stop.
%         tryagain=0;
%         stablePD=1;
%                     fprintf(' stop due to stable PD. \n ');
%      
%      end
    
%        if k>10 &&  (max(M0CVFit(k-10:k))-M0CVFit(k))./M0CVFit(k)<0.001
%             % If PD stable to within 1 percent, stop.
%         tryagain=0;
%         stableM0=1;
%                     fprintf(' stop due to stable M0. \n ');
%      
%      end
    
    % We have to stop some time if it's not convarging
    if k == maxLoops
        tryagain=0;
        if plotFlag==1
            fprintf(' rich  % 0.f  iterations. Stop before convarge \n :',maxLoops);
        end
    end
end

if plotFlag==1
    mrvNewGraphWin;
    semilogy(1:k,PDchange(1:k),'o');
    xlabel('steps'); ylabel('Change in PD')
    
    mrvNewGraphWin;
    semilogy(1:k,M0change(1:k),'o');
    xlabel('steps'); ylabel('Change in M0')
    
    mrvNewGraphWin;
    semilogy(1:k,M0Fit(1:k),'bo',1:k,M0CVFit(1:k),'ro');
    xlabel('steps'); ylabel('Change in M0 ')
    legend('Fit','CV')
end

%% Outputs

% Make a list of outputs:
OutPut.PD = PDn; % the final PD values
OutPut.Gn = Gn; % the final coils gain values
OutPut.g  = g; % the final coils gain coeficents  (G = pBasis*g)
OutPut.PDchange = PDchange;  % a vector of the PD change in each step
OutPut.M0change = M0change;  % a vector of the PD change in each step
OutPut.LastLoop = k;  % a vector of the PD change in each step
OutPut.M0Fit=M0Fit;
OutPut.M0CVFit=M0CVFit;
OutPut.NumOfIter = k;  % number of iteration

OutPut.stablePD=stablePD;
OutPut.stableM0=stableM0;
OutPut.convarge=convarge;
%keyboard
 end

