function [I] = information_matrix(X, theta, G, C, K)

[n, q] = size(X);
[~, numT ] = size(theta);
W = cell(1, numT+1);
% Pre allocate memory
for k = 1:(numT+1), W{k} = zeros(n); end

W{numT+1} = ( K\speye(n) - C*(G\C') );

for i = 1:numT
    W{i} = 1/theta(i) * ( sq_dist(X(:,i)') .* K );
    W{i} = W{i} * W{numT + 1};
end

I = zeros(numT+2);
for i = 2:(numT+2)
    for j = i:(numT+2)
        I(i,j) = sum(sum(W{i-1} .* W{j-1}'));
    end
end

I(1,1) = n-q;
for i = 2:(numT+2)
    I(1,i) = trace(W{i-1});
end

% For symmetry    
I = I + triu(I,1);
