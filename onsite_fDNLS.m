N     = 10;
eps   = 1; 
w     = 1;
alpha = 1;
J     = 1;

q = zeros(1,N);

q[1] = sqrt(w + 2*eps*J);

for n = 2:N
    temp = 0;
    for m = 1:(n-1)
        temp = temp + (JOperator(n-m) + JOperator(n+m))*q[m];
    end
    q[n] = ( eps*(JOperator(n,alpha)*q[1] ) )
end


function [coupling_weight] = JOperator(n,alpha)
    coupling_weight = 1/abs(n)^(1+alpha);
end