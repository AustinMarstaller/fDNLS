function onsite = onsite_sequence(w,eps,s,N,p)
onsite = zeros(p,1);  % Initialized initial condition
ze = zeta(s);
interval = (1:N)';
%% Construct the Onsite sequence
q0 = sqrt(w + 2*eps*ze); % Equation 3.1
q = zeros(N,1);
q(1) = (eps*q0)/(eps*(2*ze-2^(-s))+w);
    for n = 2:N
        sum1 = sum((abs(n-interval(1:(n-1))).^(-s) + abs(n+interval(1:(n-1))).^(-s)).*q(interval(1:(n-1))));
        q(n) = (eps*((q0)/(n^s) + sum1))/(eps*(2*ze-(2*n)^(-s)) + w);
    end

onsite(1:N) = q(N:(-1):1);
onsite(N+1) = q0;
onsite(N+2:p) = q(1:N); %Initial condition defined.
end