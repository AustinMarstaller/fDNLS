tic
%% Time-evolution of fDNLS
function [] = evolution_ode_repeating_alpha(slurm_id)
    %% Given the Slurm array job id, read the parameters.txt file
    % The job id serves as index for the parameters. e.g, slurm_id=2 grabs parameter set on line 2
    fileID = fopen('parameters.txt');
    content = textscan(fileID,'%.2f %.2f %.2f %.2f');
    fclose(fileID);
    
    epsVec = content{1}; % Generalized coupling power
    vVec   = content{2}; % Fractional power of Laplacian
    wVec   = content{3}; % Temporal frequency
    aVec   = content{4}; % Galilean boost
    display("Passing: eps = "+epsVec(slurm_id)+", a="+aVec(slurm_id)+", w="+wVec(slurm_id)+", v="+vVec(slurm_id)+", and ID="+slurm_id)
    % Pass eps, a, w, v to the main() function. 
    main( epsVec(slurm_id), aVec(slurm_id), wVec(slurm_id), vVec(slurm_id), slurm_id );     
end   

function []= main(eps, a, w, v, slurm_id)
%% Parameters
a = 7; % Fractional power of Laplacian
s = 1+a; 
h = 1; % Mesh parameter. Set it to one for inherently discrete systems
w = 1; % Angular frequency
v = 1; % Galilean boost
%eps = 2^a*gamma(s/2)/( sqrt(pi)*abs(gamma(-a/2)) )*h^(-a); %e.g. one can set h=1
eps = 1;
N = 50; % Grid
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
%m*z+((z.^2+(J*z).^2)).*(J*z)
F = @(t,z) eps*operator_matrix*z + ((z.^2+(J*z).^2)).*(J*z); % Dynamics
jacobian = @(t,z) eps*operator_matrix + 2*diag(z.*(J*z)) + J*diag(3*(z.*z) + (J*z).*(J*z)); % Jacobian

%% Initial condition: Not necessary to use the asymptotic solutions
initial_condition = onsite_sequence(w,eps,s,N,p);
% onsite = onsite_sequence(w,eps,s,N,p);
% initial_condition(1:p) = (0.5)^(1/4)*exp(-(0.5)*abs( (1:p)-(N+1) ));

% Y = sqrt(2*w)*sech(sqrt(w)*x); %Initial condition2

initial_drifted = [cos(v*x);sin(v*x)].*[initial_condition(1:p);initial_condition(1:p)];

%% Run ode23s
options = odeset('Jacobian', jacobian,'RelTol',1e-8,'AbsTol',1e-8);
[tSol,YSol] = ode23s(F,tspan,initial_drifted,options);
Real_YSol = YSol(:,1:p);
Imaginary_YSol = YSol(:,p+1:end);
%% Plot
close all
figure
contourf(x,tSol,Real_YSol.^2+Imaginary_YSol.^2,'EdgeColor','none')
title("Intensity","Parameters: \epsilon="+eps+", \alpha ="+a+", w="+w+", v="+v)

xlabel("Grid: N="+N), ylabel('time'), zlabel('log(|u|^2)')
colorbar
%%
figure
contourf(x,tSol,Real_YSol.^2+Imaginary_YSol.^2,[0.5,3],'EdgeColor','none')
title("Intensity","Parameters: \epsilon="+eps+", \alpha ="+a+", w="+w+", v="+v)

xlabel("Grid: N="+N), ylabel('time'), zlabel('log(|u|^2)')
colorbar

figure
contourf(x,tSol,log(Real_YSol.^2+Imaginary_YSol.^2),'EdgeColor','none')
title("log(Intensity)","Parameters: \epsilon="+eps+", \alpha ="+a+", w="+w+", v="+v)

xlabel("Grid: N="+N), ylabel('time'), zlabel('log(|u|^2)')
colorbar

%% Is L2 norm conserved? norm(YSol(i:))
clear mass
for i = 1:size(tSol), mass(i) = h*norm(YSol(i,:),2)^2; end

figure

amplitude=Real_YSol.^2+Imaginary_YSol.^2;
[argvalue, argmax] = max(amplitude');
plot(tSol,argmax-(N+1),'LineWidth',2)
title("Position of intensity"),
ylabel("ArgMax"), xlabel("time")
axis padded
%% Record time
t0 = toc;

end

