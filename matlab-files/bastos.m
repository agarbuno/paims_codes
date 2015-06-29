% Toy example from "Diagnostics for Gaussian process emulators". 

function f=bastos(x)

X1 = x(:,1);
X2 = x(:,2);

f = (1-exp(-0.5./X2)) .* (( 2300*X1.^3 + 1900*X1.^2 + 2092*X1 + 60)  ./...  
    (100*X1.^2 + 500*X1.^2 + 4*X1 + 20));