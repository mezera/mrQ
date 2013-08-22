%% Evaluate the pd solution without noise

% Try for 1d case first
x = 1:5;

% Some simple parameters
G1x = 2;
G2x = 3;
K = 5;

% These are the gains as a function of x for the two coils
r1 = (1 + x(:).*G1x) 
r2 = (K + x(:).*G2x)
% This is a convenient matrix way to calculate the r1 and r2.  Might be
% helpful some day.
% r1 = [ ones(size(x(:))) ,x(:) ]*[1, G1x]'
% r2 = [ ones(size(x(:))) ,x(:) ]*[K, G2x]'

% Write out the ratio
r = r1./r2

lhs = -1*ones(size(x(:)));
rhs = [ x(:), -1*ones(size(x(:))).*r(:) ,-1*x(:).*r(:)]

% Condition number is scary.  Not happy
svd(rhs)

% Still, we got a good solution
est = rhs\lhs

%% Try a square case
x = 1:10;

% Some simple parameters
G1x = 2; G1xx = 0.5;
G2x = 3; G2xx = 0.7;
K = 5;

% These are the gains as a function of x for the two coils
r1 = (1 + x(:).*G1x + (x(:).^2) .* G1xx);
r2 = (K + x(:).*G2x + (x(:).^2) .* G2xx);
r = r1 ./ r2;

lhs = -1*ones(size(x(:)));
rhs = [ x(:), x(:).^2, -1*ones(size(x(:))).*r(:) ,-1*x(:).*r(:), -1*(x(:).^2).*r(:)]

% Condition number is scary.  Not happy
cond(rhs)

% Still, we got a good solution
est = rhs\lhs

%% Working up to the 2D case


%
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

% Pretty close, no?
max(abs(est(:) - [G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ]'))

%% Add some noise and try again


% Seems OK to 1e-3 for this case.  But at 1e-2, things go bad.  This is a
% very high SNR level.  Uh oh.
noiseLevel = 1e-3;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;
mean(r1)/noiseLevel

r = r1Noise ./ r2Noise;

% This matrix has a low condition number.  We are going to have to fight
% this because of the problem that makes with the noise.  We may have to go
% to adding up coils (creating a virtual coil) with the svd of the data.
rhs = [X(:), X2, ...
    Y, Y2, XY, ...
    -r, -r.*X(:), -r.*X2, ...
    -r.*Y, -r.*Y2, -r.*XY];
cond(rhs)

est = rhs\lhs

% Bad noise characteristics?
max(abs(est(:) - [G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ]'))

plot(est(:),[G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ],'o')
axis square
axis equal
grid on
axis([0 6 0 6])

%% Now try some different parameters, hoping for better condition number
nSamples = 30;
[X Y] = meshgrid(1:nSamples,1:nSamples);
X = X(:);
Y = Y(:);
X2 = X(:).^2;
Y2 = Y(:).^2;
XY = X(:).*Y(:);

%
% Some simple parameters
G1x = -40; G1xx = 0.5; G1y = 2; G1yy = 1.5;
G2x = 30; G2xx = 0.7; G2y = 30; G2yy = 2.5;
G1xy = -1.2;
G2xy = -1.7;
K = 5;

r1 = 1 + X(:)*G1x + X2*G1xx + Y*G1y + Y2*G1yy + XY*G1xy;
r2 = K + X(:)*G2x + X2*G2xx + Y*G2y + Y2*G2yy + XY*G2xy;
r = r1./r2;
subplot(2,1,1), imagesc(reshape(r1,nSamples,nSamples)); colormap(gray)
subplot(2,1,2), imagesc(reshape(r2,nSamples,nSamples)); colormap(gray)

%%
lhs = -1*ones(size(X(:)));
rhs = [X(:), X2, Y, Y2, XY, -r, -r.*X(:), -r.*X2, -r.*Y, -r.*Y2, -r.*XY];

% Not a very good condition number, again.
cond(rhs)
est = rhs\lhs

% Pretty close, no?
max(abs(est(:) - [G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ]'))

%% Add some noise and try again

% Same story, though for more (x,y) samples the noise tolerance improves a
% small amount.  It is still really crummy, however.

% Seems OK to 1e-3 for this case.  But at 1e-2, things go bad.  This is a
% very high SNR level.  Uh oh.
noiseLevel = 1e-1;
r1Noise = r1 + randn(size(r1))*noiseLevel;
r2Noise = r2 + randn(size(r2))*noiseLevel;

% SNR in dB.  Very high, and yet the estimates are crummy for many of the
% gain functions.
20*log10(mean(r1)/noiseLevel)

r = r1Noise ./ r2Noise;

% This matrix has a low condition number.  We are going to have to fight
% this because of the problem that makes with the noise.  We may have to go
% to adding up coils (creating a virtual coil) with the svd of the data.
rhs = [X(:), X2, ...
    Y, Y2, XY, ...
    -r, -r.*X(:), -r.*X2, ...
    -r.*Y, -r.*Y2, -r.*XY];
cond(rhs)

est = rhs\lhs

% Bad noise characteristics?
max(abs(est(:) - [G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ]'))

plot(est(:),[G1x , G1xx , G1y, G1yy, G1xy , K, G2x, G2xx, G2y , G2yy, G2xy ],'o')
axis square
axis equal
grid on
axis([-6 6 -6 6])

%% End