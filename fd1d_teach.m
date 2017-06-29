function fd1d (M,a,b,ua,ub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       FD1D
%%       for MTH 453/553
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function shows how to set up the finite difference approximation of a 
%% model two-point boundary value problem -u''=f, on the interval [a,b], with 
%% Dirichlet boundary conditions ua (left) and ub (right), and default 
%% discretization 0....M.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Notes:
%% 1) Since MATLAB vectors aren't 0-indexed, we use 1..M+1 to label the nodes.
%% 2) The user must code the expression for f below (rhsfun).
%% 3) The user must code the expression for the exact solution below (exfun).


    %% >>> Question to MP: What do you mean by this?  It seems these are inputs ua and ub?
%% and set up /change boundary conditions 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Outputs:
%% This function will output a plot of both the exact solution and the 
%% numerical solution, along with reporting the max-norm error between the two.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% M. Peszynska, 03/2005
%% © Copyright Department of Mathematics, Oregon State University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
fprintf('%g    %g\n',dx,norm(error,inf));   %% >>> Question to MP: I don't see the two-norm here?  Am I mis-understanding?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% these functions have to be coded by user %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = rhsfun(x)  %% this is the function f from the equation -u''=f
v = pi*pi*sin(pi*x);

function v = exfun(x)   %%this is the exact solution for instructive comparison
v = sin(pi*x);
