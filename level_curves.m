function level_curves(psi, Nmin, Nmax, nlevels, meanopt, X, y)

global n;

% Gaussian process settings
meanfunc = 'h';
covfunc = 'covK';

% Prior distribution for gaussian parameters
prior = 'randn';

% Grid size
N = 70;

% Generate test points in log scale
theta1 = logspace(Nmin,Nmax,N)';
theta2 = logspace(Nmin,Nmax,N)';

[Xz, Yz] = meshgrid(theta1, theta2);
params = [ reshape(Xz,N^2,1) reshape(Yz,N^2,1) ];

L = zeros(length(params), 1);

for i=1:length(L)
    L(i) = likelihood(X, y, meanfunc, covfunc, params(i,1:2), ...
        psi, meanopt);
    if imag(L(i) )
        flag = 0;
        L(i) = Inf;
    end
    if mod(i,10000) == 0
        fprintf(1,'Iteration %12i \n', i);
    end
end

fprintf(1,'\nMinimum empiric likelihood level: %4.2f\n\n', min(L));

Lmat = reshape(L, size(Xz));
% Lmat = exp(-Lmat);

hold on
% surf(Xz,Yz,Lmat);
contour3(Xz,Yz,Lmat,nlevels);
set(gca,'yscale','log');
set(gca,'xscale','log');
xlabel('$\phi_1$', 'interpreter', 'latex')
ylabel('$\phi_2$', 'interpreter', 'latex')
title(sprintf('Number of training runs = %4i', n),  'interpreter', 'latex')
% line([10^Nmin 10^Nmax], [10^Nmin 10^Nmax], 'LineStyle', '--');
hold off