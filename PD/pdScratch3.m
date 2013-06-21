%% Working out the multiple coil case.
%
% Adding some functions to create the polynomial matrices we need.
% See the pdScratch2 for more detail.  This one is moving a bit beyond
% that.
%

%% First example, 2nd order, 2d

%T Set up parameters for three coils
G1x =  2; G1xx = 0.5; G1y = 1; G1yy = 1.5; G1xy = 1.2;  K1 = 1;
G2x = -5; G2xx = -0.7; G2y = 5; G2yy = 2.5; G2xy = 1.7;  K2 = 50;
G3x = -5; G3xx = 0.2; G3y = 5; G3yy = 0.5; G3xy = 0;    K3 = 60;

%% Initialize the matrices and responses
nSamples = 15;
[pMatrix, s] = polyCreateMatrix(nSamples,2,2);
rSize = length(s);

% r1 = 1 + X(:)*G1x + X2*G1xx + Y*G1y + Y2*G1yy + XY*G1xy;
% r2 = K + X(:)*G2x + X2*G2xx + Y*G2y + Y2*G2yy + XY*G2xy;
r1 = pMatrix*[K1,G1x,G1xx,G1y,G1yy,G1xy]';
r2 = pMatrix*[K2,G2x,G2xx,G2y,G2yy,G2xy]';
r3 = pMatrix*[K3,G3x,G3xx,G3y,G3yy,G3xy]';
%
figure;
subplot(3,1,1); imagesc(reshape(r1,rSize,rSize)); axis image
subplot(3,1,2); imagesc(reshape(r2,rSize,rSize)); axis image
subplot(3,1,3); imagesc(reshape(r3,rSize,rSize)); axis image

Z = zeros(size(pMatrix));

%% Solve for coils 1 and 2
r = r1./r2;
% lhs = -1*O;
% rhs = [M1, diag(-r(:))*M2];

rhs = [pMatrix, diag(-r(:))*pMatrix];
% [U, S, V] = svd(tmp'*tmp);
% est = U(:,end)/U(1,end);

[U, d] = eig(rhs'*rhs);
est = U(:,1)/U(1,1);
sqrt(diag(d))

% cond(rhs)
% est = rhs\lhs

% Gain parameter estimates
estMatrix = reshape(est,6,2)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2]
    

%% Add noise to the measurements (r1, r2, r3) and try again

% Seems OK to 1e-3 for this case.  But at 1e-2, things go bad.  This is a
% very high SNR level.  Uh oh.
noiseLevel = 1e-2;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;
20*log10(mean(r1)/noiseLevel)

% New r, with noisy estimates
rhs = [pMatrix, diag(-r1Noise ./ r2Noise)*pMatrix];
% [U, S, V] = svd(tmp'*tmp);
% est = U(:,end)/U(1,end);

[U, d] = eig(rhs'*rhs);
est = U(:,1)/U(1,1);
% sqrt(diag(d))

% Gain parameter estimates
estMatrix = reshape(est,6,2)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2]


%% Solve for coils 1 and 3
rhs = [pMatrix, diag(-r1./r3)*pMatrix];
[U, d] = eig(rhs'*rhs);
est = U(:,1)/U(1,1);
% sqrt(diag(d))
% Gain parameter estimates
estMatrix = reshape(est,6,2)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G3x , G3xx , G3y , G3yy , G3xy , K3]

%% Now, build up the more complex matrices.
rhs = ...
    [pMatrix, diag(-r1./r2)*pMatrix, Z; ...
     pMatrix, Z, diag(-r1./r3)*pMatrix; ...
    Z, pMatrix, diag(-r2./r3)*pMatrix];
[U, d] = eig(rhs'*rhs);
est = U(:,1)/U(1,1);

% Gain parameter estimates
estMatrix = reshape(est,6,3)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2;
    G3x , G3xx , G3y , G3yy , G3xy , K3]

%% Best so far.
%
% Notes - there is one really big value in the rhs that I don't understand.
% Adding the M2, M3 condition isn't working.
% Keep thinking.  Shouldn't adding more coils help?
noiseLevel = 1e-0;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;
r3Noise = r3 + randn(size(r3))*noiseLevel;
%
20*log10(mean(r1)/noiseLevel)

rhs = ...
    [pMatrix, diag(-r1Noise./r2Noise)*pMatrix, Z; ...
     pMatrix, Z, diag(-r1Noise./r3Noise)*pMatrix; ...
    Z, pMatrix, diag(-r2Noise./r3Noise)*pMatrix];
[U, d] = eig(rhs'*rhs);
est = U(:,1)/U(1,1);

% Gain parameter estimates
estMatrix = reshape(est,6,3)';
estMatrix = circshift(estMatrix,[0,-1]);
estMatrix
[G1x , G1xx , G1y , G1yy , G1xy,  K1;
    G2x , G2xx , G2y , G2yy , G2xy,  K2;
    G3x , G3xx , G3y , G3yy , G3xy , K3]
    
%% End