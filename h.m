% Mean function of the Gaussian process

function h = h(x,opt)

[n,m] = size(x);

switch opt
	case 1	% h(x) = 1
        h = ones(n,1);
    case 2	% h(x) = x
        h = [ones(n,1), x];
    case 3	% h(x) = 0
        h = zeros(n,1); 
    otherwise
        h = [];
        fprintf(1,'Not a valid option for h(x)\n')
end


