%Long range interaction on the lattice [-N,N] with vanishing BC on \pm N+1
%Cf Finite Difference Scheme Paper by Hong
%p = dim = 2N+1

function M = LRI_matrix(N,a)
%% Parameters
    p = 2*N+1; 
    s = 1+a; %shorthand
%% Define the diagonal part of the matrix M.
    M = 2*zeta(s)*speye(p);

%% Define the dense part
    dense = zeros(p);
    dense(1,2:p) = -(1:(p-1)).^(-s);
    for i = 2:(p-1)
        for j = (i+1):p
            dense(i,j) = dense(i-1,j-1);
        end
    end
M = M + dense + dense';
end