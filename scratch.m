%% Linear testing

x = 1:10; x = x(:);
y = 2*x;  y = y(:);

W = eye(length(x)) - (y'*y / norm(y,2));
x'*W*x
