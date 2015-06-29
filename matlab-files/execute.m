%% Driver to execute the parallel aims-opt algorithm

%% Setting the workspace and parameters of the algorithm
% clc;
clear;

rng default;
global n
n = 20;
func = 'branin';
[X,y] = startup(n, func);
load('branin_sample') % Simulator data of branin function
% load('bastos_sample') % Simulator data of bastos function
meanfunc = 'h';
covfunc = 'covK';
prior = 'randn';

gamma = .5;
alpha = 0.1;
meanopt = 2;

fprintf(1, '===================================================\n');
fprintf(1, 'Problem: ................................. %s\n', func);
fprintf(1, 'Dimension ................................ %3i\n', size(X,2));
fprintf(1, 'Training runs ............................ %3i\n', n);

%% Running aims opt
N = 2000;
lb = 10^-12;
c = 2.38/sqrt(size(X,2)+1);
fprintf(1, 'Spread parameter ......................... %1.4f\n', c);
fprintf(1, '===================================================\n\n');

% rng default

tic
    [theta, psinew, Hnew, k, w, Theta, Accep, Tvec] = parallel_aims_opt(X, y, gamma, alpha, N, c, ...
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
G = H'*C;
betahat = G\(H'*a)
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2)

%% Plotting samples 
print = false;

% psiopt = min(rescale_psi);

figure(1); clf;
level_curves(min(rescale_psi),-5,8,30,meanopt,X, y);
figure(1);
set(gca, 'FontSize', 15)
pause
hold on 
for i = 1:(size(Theta,2)/2)
    title1 = sprintf('Number of training runs = %4i', n);
    title2 = sprintf('. Temperature = %4.2f', Tvec(i));
    title(strcat(title1, title2),  'interpreter', 'latex')
    h1 = scatter(exp(Theta(:,2*i-1)), exp(Theta(:,2*i)));
    axis([10^-5 10^8 10^-5 10^8]);
    
    if print 
    % Printing options
        set(gcf, 'PaperPosition', [0.6350 6.3500 20.3200 15.2400]);
        set(gcf, 'PaperSize', [21.0000 29.7000]); 
        titlefile = sprintf('~/Dropbox/Phd/Research/Reports/Adaptive_intlik/figures/branin_wiresamp_woadap%03d.eps' ,i);
        saveas( gcf, titlefile, 'eps2c')
    end 
    pause
    if i < (size(Theta,2)/2),   delete(h1),     end;
end
scatter(exp(thetaopt(1)), exp(thetaopt(2)), 50, 'filled', 'black')
hold off