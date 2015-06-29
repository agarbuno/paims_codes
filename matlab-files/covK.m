function K = covK(X, theta, sigma, z)

% Uses the sq_dist function from Rasmussen et. al.

% Covariance function 
% Inputs: X: Matrix with data points (N x d)
%         theta: weigths of each dimension ( 1 x d )

[n,~] = size(X);
theta = sqrt(theta);

if nargin < 4
    K = sq_dist( (X/diag(theta))');
else
    K = sq_dist( (z/diag(theta))' , (X/diag(theta))');
end

K = sigma^2 * exp(- 0.5 * K);
% Force diagonal to be 1, fix numerical computation
if nargin < 4, K(1:n+1:n*n) = 1; end
% Force numerical 0 to be zeros, avoid inversion issues
K(K < sqrt(eps)) = 0;
