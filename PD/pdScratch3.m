%% Working out the multiple coil case.
%
% Adding some functions to create the polynomial matrices we need.
% See the pdScratch2 for more detail.  This one is moving a bit beyond
% that.
%

%% First example, 2nd order, 2d

% Some simple parameters
G1x = 2; G1xx = 0.5; G1y = 1; G1yy = 1.5;
G2x = 3; G2xx = 0.7; G2y = 3; G2yy = 2.5;
G1xy = 1.2;
G2xy = 1.7;
K = 5;

pMatrix = polyCreateMatrix(nSamples,2,2);
r1 = ones(nSamples*nSamples,1) + pMatrix*[G1x,G1xx,G1y,G1yy,G1xy]';
r2 = K*ones(nSamples*nSamples,1) + pMatrix*[G2x,G2xx,G2y,G2yy,G2xy]';

% r1 = 1 + X(:)*G1x + X2*G1xx + Y*G1y + Y2*G1yy + XY*G1xy;
% r2 = K + X(:)*G2x + X2*G2xx + Y*G2y + Y2*G2yy + XY*G2xy;
r = r1./r2;

lhs = -1*ones(size(pMatrix(:,1)));
pMatrix2 = [ones(nSamples*nSamples,1), pMatrix];
rhs = [pMatrix, diag(-r(:))*pMatrix2];

% rhs = [X(:), X2, Y, Y2, XY, -r, -r.*X(:), -r.*X2, -r.*Y, -r.*Y2, -r.*XY];

cond(rhs)
est = rhs\lhs

% Good fit.
max(abs(est(:) - [G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ]'))

%%
%% Add some noise and try again


% Seems OK to 1e-3 for this case.  But at 1e-2, things go bad.  This is a
% very high SNR level.  Uh oh.
noiseLevel = 1e-2;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;
20*log10(mean(r1)/noiseLevel)

% New r, with noisy estimates
r = r1Noise ./ r2Noise;
rhs = [pMatrix, diag(-r(:))*pMatrix2];

est = rhs\lhs

% Bad noise characteristics?
max(abs(est(:) - [G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ]'))

%% Now, build up the more complex matrices.

%% End