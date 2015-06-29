function [thetanew, psinew, Hnew, accep, ratenew] = delayed_parallel_annealing_level(X, y, theta, psi, ...
        N, w, c0, rate, H, Temp, meanfunc, covfunc, meanopt, lb, level, loc)

if N == 0   
    
    thetanew = []; 
    psinew = [];
    Hnew = []; 
    accep = []; 
    ratenew = [];
    
else
    
    accep = 1;  loc_accep = 1;

    % Adapting the spread
    ratenew = rate;
    c = c0 * ratenew.^(level);

    % Instead of selecting the most likely candidate to start the chain, use
    % the current loc level
    index = loc;

    % Pre-allocate memory
    thetanew = zeros(N,size(theta,2));
    psinew = zeros(N,size(psi,2));

    params = [theta, psi];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Using the covariance of the estimates %%%%
    % C = cov(params); 
    if level == Inf 
        C = 5 * speye(size(params,2));
    else
        C = weightedcov(params, w);
    end
    
    try 
        G = chol(C)';
    catch
        fprintf(1, '... Actually just use the variances ...\n');
        G = diag(std(params));
    end
    
    % G = diag(std(params));
    % G = 5 * eye(size(params,2));
    w_params = (G\params')';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [~, Ntheta] = size(theta);
    Hnew = zeros(N,1);

    % Initial state of the Markov Chain
    m = [theta(index,:), psi(index, 1)];
    counter = 0;

    while counter == 0 || (Hnew(1,1) == Inf || Hnew(1,1) == -Inf || imag(Hnew(1,1)) ~= 0 )

        newparams = m + c * ( G * randn(Ntheta + 1, 1))';

        % Extraction of the parameters
        thetanew(1,:) = newparams(1:Ntheta);
          psinew(1,:) = newparams(Ntheta+1);

        rescale_psi = (1-lb)/(1+exp(-psinew(1,1))) + lb;

        if Temp < 0
            Hnew(1,1) = likelihood_eig(X, y, meanfunc, covfunc, exp(thetanew(1,:)),...
                    rescale_psi, meanopt);
        else
            Hnew(1,1) = likelihood(X, y, meanfunc, covfunc, exp(thetanew(1,:)),...
                    rescale_psi, meanopt);
        end
        counter = counter + 1;
    end        

    for i = 1:(N-1)

        k = randsample(length(w), 1, true, w);

        % Generate local proposal
        m = [theta(k,:), psi(k,1)];
        xiparams = m + c * ( G * randn(Ntheta + 1, 1))';

        xitheta = xiparams(1:Ntheta);
          xipsi = xiparams(Ntheta + 1);
        rescale_psi = (1-lb)/(1+exp(-xipsi)) + lb;

        % This compute an approximation to the likelihood
        if Temp < 0
            Hxi = likelihood_eig(X, y, meanfunc, covfunc, exp(xitheta), ...
                rescale_psi, meanopt);
        else
            Hxi = likelihood(X, y, meanfunc, covfunc, exp(xitheta), ...
                rescale_psi, meanopt);
        end

        % This computes the local acceptance probability
        c1 = min(1,exp(-(Hxi - H(k))/Temp) );

        % Reject local proposal with proabibility 1 - rand(1,1)
        if  (Hxi == Inf || Hxi == -Inf || imag(Hxi) ~= 0 || isnan(Hxi) )  || rand(1,1) >= c1
            thetanew(i+1, :) = thetanew(i, :);
              psinew(i+1, 1) = psinew(i, 1);

            Hnew(i+1,1) = Hnew(i,1);
        else

            loc_accep = loc_accep + 1;

            % cand means the candidate for the markov chain
            % curr means the current step in the markov chain

            c2 = min(1, exp(-(Hnew(i,1) - H)/Temp )) ;

            % Current state of the Markov Chain
            newparams = [thetanew(i,:), psinew(i,1) ];

            q_curr = 1/sqrt(c^2 * prod(diag(G).^2)) * ...
                (2*pi)^(-size(params,2)/2) * exp(- (0.5/c^2) * ...
                sum(bsxfun( @minus, w_params, (G\newparams')').^2, 2));
            phat_curr = (w .* q_curr)'* c2;

            c2 = min(1, exp(-(Hxi - H)/Temp));

            q_cand = 1/sqrt(c^2 * prod(diag(G).^2)) * ...
                (2*pi)^(-size(params,2)/2) * exp(- (0.5/c^2) * ...
                sum(bsxfun( @minus, w_params, (G\xiparams')').^2, 2));
            phat_cand = (w .* q_cand)'* c2;

            % Probability of accepting the candidate as the next step in the
            % chain
            frac2 = min(1, exp(-(Hxi - Hnew(i,1))/Temp) * ...
                phat_curr/phat_cand);

            if rand(1,1) < frac2
                thetanew(i+1, :) = xitheta;
                  psinew(i+1, 1) = xipsi;

                Hnew(i+1,1) = Hxi;
                accep = accep + 1;
            else
                % Globally rejected sample time to perform a delayed rejection
                % scheme

                % Case 1 - Random walk from the current candidate
                crumb = newparams + c0 * ( G * randn(Ntheta + 1, 1))';
                aux = 1;

                crumbtheta = xiparams(1:Ntheta);
                  crumbpsi = xiparams(Ntheta + 1);
                rescale_crumbpsi = (1-lb)/(1+exp(-xipsi)) + lb;

                Hcrumb = likelihood(X, y, meanfunc, covfunc, exp(crumbtheta), ...
                    rescale_crumbpsi, meanopt);

                c2 = min(1, exp(-(Hcrumb-H)/Temp));
                q_crumb = 1/sqrt(c0^2 * prod(diag(G).^2)) * ...
                    (2*pi)^(-size(params,2)/2) * exp(- (0.5/c0^2) * ...
                    sum(bsxfun( @minus, w_params, (G\crumb')').^2, 2));
                phat_crumb = (w .* q_crumb)'* c2;

                % Probability of accepting the candidate by coming from the
                % crumb
                frac_crumbcand = min(1, exp(-(Hxi - Hcrumb)/Temp) * ...
                    phat_crumb / phat_cand);

                frac_delayed = min(1, exp(-(Hcrumb - Hnew(i,1) )/Temp) * ...
                    aux * (1- frac_crumbcand)/(1- frac2));

                if rand(1,1) < frac_delayed
                    thetanew(i+1, :) = crumbtheta;
                      psinew(i+1, 1) = crumbpsi;

                    Hnew(i+1,1) = Hcrumb;
                    accep = accep + 1;
                else
                    thetanew(i+1, :) = thetanew(i, :);
                      psinew(i+1, 1) = psinew(i, 1);

                    Hnew(i+1,1) = Hnew(i,1);
                end

            end

        end

    end

end