%% Driver to execute the parallel aims-opt algorithm for the Nilson-Kuusk algorithm


%% Setting the workspace
% clc;
clear;

rng default;
global n
infile = 'training.txt';
outfile = 'out_5_100.txt';
X = dlmread(infile);
min_X = min(X);
max_X = max(X);
X = bsxfun(@rdivide, bsxfun(@minus, X, min(X)), (max(X) - min(X)));
y = dlmread(outfile);

meanfunc = 'h';
covfunc = 'covK';
prior = 'randn';

gamma = .5;
alpha = 0.1;
meanopt = 2;

%% Running aims opt
N = 5000;
lb = 10^-12;
c = 2.38/sqrt(size(X,2)+1);

fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', 'kuusk model');
fprintf(1, 'Dimension ................................ %3i\n', size(X,2));
fprintf(1, 'Training runs ............................ %3i\n', size(X,1));
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default;

tic
    [theta, psinew, Hnew, k, w, Theta, Accep] = parallel_aims_opt(X, y, gamma, alpha, N, c, ...
        meanfunc, covfunc, meanopt, prior, lb);
toc

fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));


%% Nugget term 
rescale_psi = (1-lb)./(1+exp(-psinew)) + lb;
psiopt = rescale_psi(Hnew == min(Hnew));
psiopt = unique(psiopt);

thetaopt = theta(Hnew == min(Hnew),:);
thetaopt = unique(thetaopt,'rows');

fprintf(1, 'Optimal nugget ............... %8.2e\n', psiopt);
if size(thetaopt,2) == 2
    fprintf(1, 'Optimal theta: [ %4.8f, %4.8f ] \n', exp(thetaopt(1)), exp(thetaopt(2)) );
else
    fprintf(1, 'Optimal theta: [ %4.8f, %4.8f, %4.8f, %4.8f, %4.8f ] \n', ...
        exp(thetaopt(1)/2), exp(thetaopt(2)/2), exp(thetaopt(3)/2), exp(thetaopt(4)/2), exp(thetaopt(5)/2) );
end

%% Computes estimations for sigma
params = [exp(theta)];
L = -Hnew;
n = size(X,1);

maxp = params(L == max(L),:);
Lmax = max(L);
maxp = maxp(1,:);

thetahat = maxp;

H = feval(meanfunc, X, meanopt);
K = feval(covfunc, X, thetahat, 1);
K = K + psiopt * eye(n);

a = K\y;
C = K\H;
betahat = (H'*C)\(H'*a)
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2)
