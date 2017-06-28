function fd1d (M,a,b,ua,ub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       FD1D
%%       for MTH 453/553
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function shows how to set up 
%% the finite difference approximation of a model 
%% two-point boundary value problem -u''=f
%% user must code expression for f (rhsfun)
%% and set up /change boundary conditions 
%% default discretization is 0....M, with M a parameter
%% but since MATLAB can't number vectors from 0, we use 1..M+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M. Peszynska, 03/2005
%% © Copyright Department of Mathematics, Oregon State University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1, 
    M = 10; 
    a = 0; b = 1;
    ua = 0;
    ub = 0;
end

dx = (b-a)/M;           %% use uniform grid
x = (a:dx:b)';

%% set up the matrix and the right hand side vector
f = rhsfun(x);

A = zeros(M+1,M+1);     %% or use sparse(M+1,M+1)
for i=2:M
    %% set up i'th row of the matrix
    A(i,i) = 2;
    A(i,i-1) = -1; 
    A(i,i+1) = -1; 
end
A = A / (dx*dx);

%% apply the Dirichlet boundary conditions to the matrix 
%% and to the right hand side
A (1,1) = 1;
A (M+1,M+1) = 1; 
f (1,1) = ua;
f (M+1,1) = ub;

%% solve the linear system;
u = A \ f;

%% plot both the analytical and the numerical solutions
exvec = exfun(x);
plot(x,exvec,'k',x,u','r*');

%% compute the error vector
error = abs(exvec - u);
%% report the error in max norm and in two discrete norm 
fprintf('%g    %g\n',dx,norm(error,inf));

%%%%%%%%%%%%%%%%%%%%%%%% these functions have to be coded by user
function v = rhsfun(x)
v = pi*pi*sin(pi*x);

function v = exfun(x)
v = sin(pi*x);

        
    
