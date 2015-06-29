%% Function to prepare the workspace for two dimensional examples.
% Input: 	n - number of training points
% 			func - example simulator

% Output: 	X - matrix with design points
%			y - outputs of the simulator

function [X,y] = create_data(n, func)
% clear;

% Return to original matlab format
format

%% To recover points from a latin hypercube
k = 2;
X = lhsdesign(n,k);

% Compute observed data
for i = 1:n
    y(i,1) = feval(func, (X(i,:)));
end