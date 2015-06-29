function L = likelihood(data, y, meanfunc, covfunc, theta, psi, meanopt)

% For numerical stability it should be better to return the negative of the 
% log likelihood

[n, q] = size(data);

H = feval(meanfunc, data, meanopt);
K = feval(covfunc, data, theta, 1);
K = K + psi * speye(n);

a = K\y;
C = K\H;
G = H'*C;

% Check this line
sigma = y'*(a - C*(G\(H'*a)))/(n-q-2);
detK = det(K);
% if detK < eps
%     detK = 0;
% end
L = (n-q)/2 * log(sigma) + 0.5 * log(detK) + 0.5 * log(det(G)); 

% This line computes the log-normal prior distribution
% L = L - log(prior_lognormal(theta, 3));

% Partition the components in likelihood to look at them individually
% L = (n-q)/2 * log(sigma);
% L =  0.5 * log(detK);
% L =  0.5 * log(det(H'*C));

% This line computes the reference prior distribution
[I] = information_matrix(data, theta, G, C, K);
L = L - log(sqrt(det(I)));

% This line adds a number to guarantee the nonnegative restriction of the log-posterior function
L = L + 400;