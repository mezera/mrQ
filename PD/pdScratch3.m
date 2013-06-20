%% Working out the multiple coil case.
%
% Adding some functions to create the polynomial matrices we need.
% See the pdScratch2 for more detail.  This one is moving a bit beyond
% that.
%

%% First example, 2nd order, 2d

%T Set up parameters for three coils
G1x = 2; G1xx = 0.5; G1y = 1; G1yy = 1.5; G1xy = 1.2;  K1 = 1;
G2x = 3; G2xx = 0.7; G2y = 3; G2yy = 2.5; G2xy = 1.7;  K2 = 5;
G3x = -5; G3xx = 0.2; G3y = 5; G3yy = 0.5; G3xy = 0;   K3 = 10;

%% Initialize the matrices and responses
nSamples = 30;
[pMatrix, O] = polyCreateMatrix(nSamples,2,2);

% r1 = 1 + X(:)*G1x + X2*G1xx + Y*G1y + Y2*G1yy + XY*G1xy;
% r2 = K + X(:)*G2x + X2*G2xx + Y*G2y + Y2*G2yy + XY*G2xy;
r1 = K1*O + pMatrix*[G1x,G1xx,G1y,G1yy,G1xy]';
r2 = K2*O + pMatrix*[G2x,G2xx,G2y,G2yy,G2xy]';
r3 = K3*O + pMatrix*[G3x,G3xx,G3y,G3yy,G3xy]';

M1 = pMatrix;
M2 = [O, pMatrix];
M3 = M2;

%% Solve for coils 1 and 2
r = r1./r2;
lhs = -1*O;
rhs = [M1, diag(-r(:))*M2];

% rhs = [X(:), X2, Y, Y2, XY, -r, -r.*X(:), -r.*X2, -r.*Y, -r.*Y2, -r.*XY];

% cond(rhs)
est = rhs\lhs

% Gain parameter estimates
estMatrix = reshape([1;est],6,2)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2]
    

%% Add some noise and try again


% Seems OK to 1e-3 for this case.  But at 1e-2, things go bad.  This is a
% very high SNR level.  Uh oh.
noiseLevel = 1e-3;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;
20*log10(mean(r1)/noiseLevel)

% New r, with noisy estimates
r = r1Noise ./ r2Noise;
rhs = [M1, diag(-r(:))*M2];

est = rhs\lhs

% Gain parameter estimates
estMatrix = reshape([1;est],6,2)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2]


%% Solve for coils 1 and 3
r = r1./r3;

lhs = -1*O;
rhs = [M1, diag(-r(:))*M3];
est = rhs\lhs

% Gain parameter estimates
estMatrix = reshape([1;est],6,2)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G3x , G3xx , G3y , G3yy , G3xy , K3]

%% Now, build up the more complex matrices.
Z = zeros(size(M2));
rhs = ...
    [M1, diag(-r1./r2)*M2, Z; ...
     M1, Z, diag(-r1./r3)*M3];
%  ...
%      Z(:,1:(end-1)), M2, diag(-r2./r3)*M3];
%  
NegOne = vertcat(lhs,lhs);
est = rhs\NegOne;

% Gain parameter estimates
estMatrix = reshape([1;est],6,3)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2;
    G3x , G3xx , G3y , G3yy , G3xy , K3]

%% So far, doesn't help.
% Keep thinking.  Shouldn't adding more coils help?
noiseLevel = 1e-2;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;
r3Noise = r3 + randn(size(r3))*noiseLevel;

Z = zeros(size(M2));
rhs = ...
    [M1, diag(-r1Noise./r2Noise)*M2, Z; ...
     M1, Z, diag(-r1Noise./r3Noise)*M3];
%  ...
%      Z(:,1:(end-1)), M2, diag(-r2./r3)*M3];
%  
NegOne = vertcat(lhs,lhs);
est = rhs\NegOne;
size(rhs'*rhs)

% Gain parameter estimates
estMatrix = reshape([1;est],6,3)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2;
    G3x , G3xx , G3y , G3yy , G3xy , K3]
    
%% End