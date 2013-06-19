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

