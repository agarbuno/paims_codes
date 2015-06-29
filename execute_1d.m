%% Driver to execute the parallel aims-opt algorithm in one dimensional examples

%% Setting the workspace
% clc;
clear;

% rng default;
global n
n = 8;
func = 'branin';

X = linspace(-4.35,4.35,n)';
% X = -5 + 10 * lhsdesign(n,1);
% Calculate observed data
y = oakley(X); 

meanfunc = 'h';
covfunc = 'covK';
prior = 'randn';

gamma = .5;
alpha = 0.05;
N = 3000;
c = 1.2;
meanopt = 2;

%% Running aims opt
lb = 10^-9;

tic
    [theta, psinew, Hnew, k, w, Theta] = parallel_aims_opt(X, y, gamma, alpha, N, c, ...
        meanfunc, covfunc, meanopt, prior, lb);
toc

fprintf(1, '\nInterval: [ %4.8f, %4.8f ] \n', min(Hnew), max(Hnew));


%% Nugget term 
rescale_psi = (1-lb)./(1+exp(-psinew)) + lb;
psiopt = rescale_psi(Hnew == min(Hnew));
psiopt = unique(psiopt);

thetaopt = theta(Hnew == max(Hnew),:);
thetaopt = unique(thetaopt,'rows');

fprintf(1, 'Optimal nugget ............... %8.2e\n', psiopt);
if size(thetaopt) > 1
    fprintf(1, 'Optimal theta: [ %4.8f, %4.8f ] \n', exp(thetaopt(1)), exp(thetaopt(2)) );
else
    fprintf(1, 'Optimal theta: [ %4.8f ] \n', exp(thetaopt(1)) );
end

%% Plot for one dimensional example

params = [exp(theta)];
L = -Hnew;

maxp = params(L == max(L),:);
% maxp = params(rescale_psi == min(rescale_psi),:);
% psiopt = min(rescale_psi)
Lmax = max(L);
maxp = maxp(1,:);

size(maxp)

thetahat = maxp(1,:);

H = feval(meanfunc, X, meanopt);
K = feval(covfunc, X, thetahat, 1);
K = K + psiopt * eye(n);

a = K\y;
C = K\H;
betahat = (H'*C)\(H'*a)
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2)

z = linspace(-5, 5, 100)';
k = feval(covfunc, X, thetahat, 1, z);
h = feval(meanfunc, z, meanopt);

yz = h*betahat + k * (K\( y - m ) );
kii = diag(feval(covfunc, z, thetahat, 1));

% c* for the covariance estimator
s = kii - diag(k * (K\k'));
% c** for the covariance estimator
dev = h - k*C;
s = s + diag(dev*((H'*C)\dev'));
s = sigmahat^2 * s;

f = [yz+2*sqrt(s); flipdim(yz-2*sqrt(s),1)];
% % Plot the function and the training rounds
figure(2); clf; 
Yz = oakley(z);
hold on; plot(z, Yz)
fill([z; flipdim(z,1)], f, [0.8 0.2 .2],'FaceAlpha',0.1,'EdgeAlpha',0.3); plot(z, yz, '.');
plot(X, y, 'rx'); 
hold off;

%%  All modes
num = length(rescale_psi);
z = linspace(-5, 5, 100)';
Y = [];

for i = 1:num
    
    if Hnew(i) <= min(Hnew) + std(Hnew)
    
        H = feval(meanfunc, X, meanopt);
        K = feval(covfunc, X, exp(theta(i,:)), 1);
        K = K + rescale_psi(i) * eye(n);

        a = K\y;
        C = K\H;
        betahat = (H'*C)\(H'*a);
        m = H*betahat;

        z = linspace(-5, 5, 100)';
        k = feval(covfunc, X, thetahat, 1, z);
        h = feval(meanfunc, z, meanopt);

        yz = h*betahat + k * (K\( y - m ) );

        Y = [Y, yz];
        
    end
    
end

figure(1); clf; 
yz = oakley(z);
hold on; 
plot(z, Y);
plot(X, y, 'rx'); plot(z, yz, 'k.', 'LineWidth', 6);
hold off;
