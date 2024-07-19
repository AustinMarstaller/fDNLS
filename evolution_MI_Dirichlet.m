function [] = evolution_MI_Dirichlet()
    % Pass a to main() function.
    close all
    
 
    for a = [1,20]
        for A = [1]
            figure 
            [tSol,YSol,p,N,x]  = main(a,A);     
        
            Real_YSol = YSol(:,1:p);
            Imaginary_YSol = YSol(:,p+1:end);
            pcolor(x,tSol,Real_YSol.^2+Imaginary_YSol.^2)
            shading interp
            colorbar, cmocean('balance')
            xlabel("Space")
            ylabel("Time")
            title("Intensity: \alpha = "+a)
            subtitle("A = "+A)
        end    
    end

end

function [tSol,YSol,p,N,x,A] = main(a,A)
    %% Solve the fDNLS with Dirichlet boundary conditions. Time-evolution of fDNLS
    % Problem is split into Re and Im parts. 

tic
warning('off','all')

%% Fixed parameters
eps=1;
w = 1;
v = 1;
s = 1+a; 
h = 1;   
N = 100; % Grid
p = 2*N+1; % Total amount of grid points

tspan = [0 50]; % Time domain
x = ((-N*h):(h):(N*h))'; % Space domain

%% Real/Im matrix
operator_matrix = sparse(2*p,2*p);
operator_matrix(1:p,(p+1):(2*p)) = LRI_matrix(N,a);
operator_matrix = operator_matrix - operator_matrix'; %Riesz frac Laplacian discretized (zero exterior BC)

% Symplectic matrix
J = sparse(2*p,2*p);
J(1:p,p+1:2*p) = -speye(p);
J = J - J';

%% ODE Problem defined
F = @(t,z) eps*operator_matrix*z + ((z.^2+(J*z).^2)).*(J*z); % Dynamics
jacobian = @(t,z) eps*operator_matrix + 2*diag(z.*(J*z)) + J*diag(3*(z.*z) + (J*z).*(J*z)); % Jacobian

%% Initial condition
eps=1;
%initial_condition = zeros(2*p,1);
%initial_condition(100) = A;
onsite = onsite_sequence(w,eps,s,N,p); % Execute the onsite_sequence script to generate the Onsite asymptotic sequence
initial_condition(1:p) = onsite;
%for i = 1:2*p
%initial_condition(i) = 1;
%end
% Y = sqrt(2*w)*sech(sqrt(w)*x); %Initial condition2

%initial_drifted = [cos(v*x);sin(v*x)].*[onsite;onsite]; % Galilean drift 
initial_drifted = [cos(v*x);sin(v*x)].*[initial_condition;initial_condition]; % Galilean drift 

%% Run ode23s
options = odeset('Jacobian', jacobian,'RelTol',1e-3,'AbsTol',1e-6); % Pass the Jacobian. Significant speed increase using this.
[tSol,YSol] = ode23s(F,tspan,initial_drifted,options); % Solve the problem
Real_YSol = YSol(:,1:p);
Imaginary_YSol = YSol(:,p+1:end);
%% Record time
toc
end


