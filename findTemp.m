function [tempnew, flag, counter] = findTemp(H, Temp, gamma, N)

% This function helps calculate the next value for beta
% the parameter of annealing in the next level.

TOL = 10d-8;
beta = 1/Temp;
a = beta;
b = 100000000;
flag = 0;

betanew = (a + b)/2;
k = 1/(N*gamma);

fa = (sum(exp(- 2 * H * 0)) / sum(exp(- H * 0))^2) - k;
fc = (sum(exp(- 2 * H * (betanew - beta) )) / sum(exp(- H * (betanew - beta) ))^2) - k;

Nmax = 100;
counter = 0;

fprintf(1, '\nBisection iterations ..................... ')

while counter <= Nmax
    betanew = (a + b)/2;
    fc =  (sum(exp(- 2 * H * (betanew - beta) )) / sum(exp(- H * (betanew - beta) ))^2) - k;
    
    if fc == 0 || (b - a)/2 < TOL
        flag = 1;
        break
    end
    counter = counter + 1;
    
    if sign(fc) == sign(fa) 
        a = betanew;
    else
        b = betanew;
    end
    
    fa = (sum(exp(- 2 * H * (a - beta) )) / sum(exp(- H * (a - beta) ))^2) - k; 
    
end

tempnew = 1/betanew;

fprintf(1, 'Done \n')
fprintf(1, '\tTemperature ...................... %4.3f\n', 1/betanew)

