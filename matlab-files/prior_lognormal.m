function y = prior_lognormal(x, sigma)

% This function calculate the pdf for iid lognormal random variables in a
% row vector. If feed with a matrix it returns a vector with the pdf for
% every row in matrix x.

y = prod((1./(x*sigma*sqrt(2*pi))) .* exp(- (0.5/sigma^2) * log(x).^2));