%% Draw test dataset

infile = '~/Dropbox/Datafun/model2/in_5_150.txt';
outfile = '~/Dropbox/Datafun/model2/out_5_150.txt';
Xtry = dlmread(infile);
Xtry = bsxfun(@rdivide, bsxfun(@minus, Xtry, min_X), (max_X - min_X));
Ytry = dlmread(outfile);

ktry = size(Xtry,1);

%% Computes estimations for sigma
params = exp(theta);
L = -Hnew;
n = size(X,1);

maxp = params(L == max(L),:);
Lmax = max(L);
maxp = maxp(1,:);

thetahat = maxp;

H = feval(meanfunc, X, meanopt);
K = feval(covfunc, X, thetahat, 1);
K = K + psiopt * eye(n);

a = K\y;    C = K\H;    G = H'*C;
betahat = G\(H'*a);
m = H*betahat;

sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2);

%% Plot standarized residuals

k = feval(covfunc, X, thetahat, 1, Xtry);
h = feval(meanfunc, Xtry, meanopt);

% Estimated mean
mean_y = h*betahat + k * (K\( y - m ) );
RMSE = sqrt(mean((mean_y - Ytry).^2))


% Error prediction
Kchol = chol(K)'; Kk = (Kchol\k')';
c_star = ones(size(k,1),1) - sum(Kk .* conj(Kk),2);
A = (h - k * (K\H)); Gchol = chol(G)'; GA = (Gchol\A')';
c_star = c_star - sum(GA .* conj(GA),2);
var_y = sigmahat * c_star;

% Standarized independent residuals
res_ind = (Ytry - mean_y)./sqrt(var_y);

% Plotting the residuals
% Plotting the residuals
figure(1); clf;
ylabel('Standarized Residuals', 'interpreter', 'latex')
subplot(2,1,1);
hold on
plot(res_ind, '.', 'MarkerSize',10)
title('Residual plot using MAP' ,  'interpreter', 'latex')
plot(get(gca,'xlim'), [3 3], '--k'); 
plot(get(gca,'xlim'), [-3 -3], '--k'); 
set(gca, 'FontSize', 15)
axis([0 Inf -10 10])
hold off

%plot(mean_y, res_ind, '.')
%plot(Xtry(:,1), res_ind,'.')
%plot(Xtry(:,2), res_ind,'.')

%% Using the mixture

y_mix = zeros(size(Xtry,1),1);
w_new = 1/size(theta,1) * ones(size(theta,1),1);
meanyz = zeros(size(Xtry,1), size(theta,1));
varyz = zeros(size(Xtry,1), size(theta,1));

% This for iterates through all the samples from paims
for i = 1:size(theta,1)
    K = feval(covfunc, X, exp(theta(i,:)), 1);
    K = K + rescale_psi(i,:) * speye(n);

    a = K\y;    C = K\H;    G = H'*C;    betahat = G\(H'*a);
    m = H*betahat;

    sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2); 
    
    k = feval(covfunc, X, exp(theta(i,:)), 1, Xtry);
    yz = h*betahat + k * (K\( y - m ) );
    
    y_mix = y_mix + w_new(i) * yz;
    meanyz(:,i) = yz;
    
    Kchol = chol(K)'; Kk = (Kchol\k')';

    % This is an efficient way to compute the diagonal of a product of
    % matrices. In comments are the past computations, more readable but
    % way more expensive.
    %     c_star = ones(size(k,1),1) - diag(k*(K\k'));
    c_star = ones(size(k,1),1) - sum(Kk .* conj(Kk),2);

    A = (h - k * (K\H)); Gchol = chol(G)'; GA = (Gchol\A')';
    c_star = c_star - sum(GA .* conj(GA),2);
       
    varyz(:,i) = sigmahat * c_star;
end

RMSE = sqrt(mean((y_mix - Ytry).^2))
var_mix = bsxfun(@plus, bsxfun(@minus, meanyz, y_mix).^2, varyz) * w_new;

% Standarized independent residuals
res_ind_mix = (Ytry - y_mix)./sqrt(var_mix);

% Plotting the residuals
subplot(2,1,2);
hold on
plot(res_ind_mix, '.', 'MarkerSize',10)
xlabel('Index', 'interpreter', 'latex')
title('Residual plot using mixture model' ,  'interpreter', 'latex')
plot(get(gca,'xlim'), [3 3], '--k'); 
plot(get(gca,'xlim'), [-3 -3], '--k'); 
set(gca, 'FontSize', 15)
axis([0 Inf -10 10])
hold off


%plot(mean_y, res_ind_mix, '.')
%plot(Xtry(:,1), res_ind_mix,'.')
%plot(Xtry(:,2), res_ind_mix,'.')

%% Comparing residuals

% figure(3); clf;
% plot(res_ind, res_ind_mix, '.')
% xlabel('Standarized Residuals MLE', 'interpreter', 'latex')
% ylabel('Standarized Residuals mixture model', 'interpreter', 'latex')
% set(gca, 'FontSize', 15)
% 
% %% Obtain the correlation matrix
% C_matrix = zeros(size(Xtry,1));
% 
% for i = 1:size(theta,1)
%     K = feval(covfunc, X, exp(theta(i,:)), 1);
%     K = K + rescale_psi(i,:) * speye(n);
% 
%     a = K\y;    C = K\H;    G = H'*C;
%     
%     sigmahat = y'*(a - C*((H'*C)\(H'*a)))/(n-size(X,2)-2); 
%     
%     k = feval(covfunc, X, exp(theta(i,:)), 1, Xtry);
%     
%     Kchol = chol(K)'; Kk = (Kchol\k')';
%     c_star = feval(covfunc, Xtry, exp(theta(i,:)), 1) - Kk*Kk';
%    
%     A = (h - k * (K\H)); Gchol = chol(G)'; GA = (Gchol\A')';
%     c_star = c_star - GA*GA';
%        
%     C_matrix = C_matrix + w_new(i) * sigmahat * c_star;
% end
% 
% aux = bsxfun(@minus, meanyz, y_mix);
% C_matrix = aux * (diag(w_new) * aux') + C_matrix;
% C_matrix = (C_matrix + C_matrix')/2;
% 
% %% Plot pivoted cholesky errors
% 
% [Chol, P] = cholp(C_matrix);
% 
% res_chol_mix = (Chol*P')\(y_mix - Ytry);
% % Plotting the residuals
% figure(4), clf;
% plot(1:150, res_chol_mix, '.')
% % plot(P*(1:150)', res_chol_mix, '.')
% xlabel('Index', 'interpreter', 'latex')
% ylabel('Standarized Residuals', 'interpreter', 'latex')
% title('Cholesky errors using mixture model' ,  'interpreter', 'latex')
% set(gca, 'FontSize', 15)

