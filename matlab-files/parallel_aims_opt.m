function [theta, psi, H, j, w, Theta, Accep, Tvec] = parallel_aims_opt(X, y, gamma, alpha, N,...
    c0, meanfunc, covfunc, meanopt, prior, lb)

[~,d] = size(X);

% Rate of reduction for spread, learned online through the algorithm.
% rate = 1 * rand(N, 1);
rate = .5 * ones(1,1);

%   Annealing level j
j = 0; 
Temp = Inf; Tvec = [];
delta = 0;
maxSamp = 20;
flag = false;

% Pior chain
mag = .75;
counter = 0;
H = zeros(N(1),1);

% Check if we have a previous sample to start from
while (delta == 0 || isnan(delta)) && counter < maxSamp
    if counter > 0, 
        fprintf(1, 'I am really trying ....... %3i', counter); 
        min(H)
    end
    
    if strcmp(covfunc, 'covKiso')
        theta = mag * feval(prior,N(1),1);
    else
        theta = -7 + (14) * rand(N, d);
    end

    psi = 10 * randn(N(1), 1);

    H = zeros(N(1),1);
    parfor i = 1:N(1)
        rescale_psi = (1-lb)/(1+exp(-psi(i,1))) + lb;
        H(i) = likelihood(X, y, meanfunc, 'covK', exp(theta(i,:)),...
            rescale_psi, meanopt );
    end
    
    indexinf = (H == Inf); indexminf = (H == -Inf); indeximag = (imag(H)~=0);
    psi = psi(~(indexinf+ indexminf+ indeximag),:);
    theta = theta(~(indexinf+ indexminf+ indeximag),:);
    H = H(~(indexinf+ indexminf+ indeximag),:);
    
    delta = sqrt(mean((H - mean(H)).^2)) / mean(H);
    counter = counter + 1;

end

Accep = []; Theta = [];

covfunc = 'covK';
fprintf(1, 'Initial samples .......................... %4i', length(H));
if N > length(H)
    indexes = randsample(length(H), N - length(H));
    H = [H ; H(indexes,:)];
    theta = [theta ; theta(indexes,:)];
    psi = [psi ; psi(indexes,:)];
end

deltatarget = alpha * delta;
fprintf(1, '\nTarget delta ............................. %4.4e\n', deltatarget);

% Pre allocating memory
Ntheta = d;
thetanew = zeros(N, Ntheta);
psinew = zeros(N,1); accep = 0;
ratenew = zeros(N,1); Hnew = zeros(N,1);

while ~flag && (delta > deltatarget && j < 20 )
    
    fprintf(1, '\nStarting annealing level: %4i', j+1);
    Tempnew = findTemp(H, Temp, gamma, N);
    
    % This ends the algorithm with a sample of the posterior distribution
%     if Tempnew < 1
%         flag = true;
%         Tempnew = 1;
%     end    
    
    w = exp ( - H * (1/Tempnew - 1/Temp) );
    w = w / sum(w);

    % Generate samples from wn to find how many times each chain will be updated. 
    Ns = zeros(1,N);    
    cx = [0; cumsum(w)];
    for i = 1:N
        n = 1;
        u = rand;
        while u > cx(n)
            n = n + 1;
        end
        Ns(n-1) = Ns(n-1) + 1;
    end  
    
    fprintf(1, 'Coefficient of Correlation ............... %3.4f\n', 100*std(w)/mean(w));
    fprintf(1, 'Effective Sample Size .................... %4.4e\n', 1/sum(w.^2));
    fprintf(1, 'Temperature difference ................... %4.2f\n', (Temp - Tempnew)/Temp);
    fprintf(1, 'Entering level ...........................%4i\n', j+1);
    
    % Discard the chains not worth growing
    theta = theta(Ns>0,:); psi = psi(Ns>0,:);
    H = H(Ns>0,:); w = w(Ns>0,:); Ns = Ns(Ns>0);
    
    % Run the chain in parallel
    parfor i = 1:length(Ns)
        [ thetanew_sampled{i}, psinew_sampled{i}, Hnew_sampled{i}, ...
            accep_sampled{i}, ratenew_sampled{i} ] = ...
            delayed_parallel_annealing_level(X, y, theta, psi, ...
            Ns(i), w, c0, rate, H, Tempnew, meanfunc, covfunc, meanopt, lb, j, i);
    end
    
    % Gather all the new samples together    
    n1 = 1;
    for i = 1:length(Ns)
        if isempty(thetanew_sampled{i}) == 1
        else
            n2 = n1 + length(thetanew_sampled{i}(:,1));
            thetanew(n1:n2-1,1:Ntheta) = thetanew_sampled{i};
            psinew(n1:n2-1,:) = psinew_sampled{i};
            Hnew(n1:n2-1,:) = Hnew_sampled{i};
            accep = accep + accep_sampled{i};
            ratenew(n1:n2-1,:) = ratenew_sampled{i};
            n1 = n2;
        end
    end
    
    ratenew = mean(ratenew); accep = accep/N;
    
    deltanew = sqrt(mean((Hnew - mean(Hnew)).^2)) / mean(Hnew);
    
    theta = thetanew; psi = psinew;
    deltachange = abs((deltanew - delta)/delta);
    H = Hnew; delta = deltanew; j = j + 1;
    
    % Resampling idea, not used in the final version of the algorithm
    % [theta, psi, H, w, ratenew] = resampling(H, Temp, Tempnew, N, gamma, theta, psi, ...
    %     X, y, meanfunc, covfunc, c0, ratenew, j-1, lb, meanopt, accep);
    
    rate = ratenew;
        
    fprintf(1, 'Global Acceptance rate ................... %5.2f \n', accep * 100);
    fprintf(1, 'Rate of reduction ........................ %5.5f \n', ratenew);
    fprintf(1, 'Delta value .............................. %4.4e\n', delta);
    fprintf(1, 'Delta change ............................. %4.4e\n', deltachange);
    
    Accep = [Accep, accep];
    Theta = [Theta, theta];
    
    Temp = Tempnew; Tvec = [Tvec; Temp];

end    

w = exp ( - H * (1/Temp) );
w = w / sum(w);
fprintf(1, '\nLast Effective Sample Size ............... %4.4e\n', 1/sum(w.^2));


%  This might be needed if optimisation is not pursued.

% if flag 
%     w = exp ( - H * (1/Tempnew - 1/Temp) );
%     w = w / sum(w);
%     [ thetanew, psinew, Hnew, accep, ratenew] = annealing_level(X,y, theta, psi, ...
%         N, w, c0, rate, H, 1, meanfunc, covfunc, meanopt, lb, j);
%     theta = thetanew; psi = psinew;
%     
%     deltachange = abs((deltanew - delta)/delta);
%     H = Hnew; delta = deltanew; j = j + 1;
%     
%     [theta, psi, H, w, ratenew] = resampling(H, Temp, Tempnew, N, gamma, theta, psi, ...
%         X, y, meanfunc, covfunc, c0, ratenew, j-1, lb, meanopt, accep);
%     
%     Accep = [Accep, accep];
%     Theta = [Theta, theta];
%     
%     Temp = 1; Tvec = [Tvec; Temp];
% end
