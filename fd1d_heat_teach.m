function [maxerr] = fd1d_heat (M,a,b,ua,ub,T,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finite difference approximation of a simple 
%% parabolic two-point boundary value problem u_t-u_{xx} = f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% M. Peszynska,  2005/4
%% Copyright Department of Mathematics, Oregon State University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%% INSTRUCTIONS: the code can be executed as script or retyped
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% user must code expression for f (rhsfun)
%% for initial condition initfun()
%% and set up /change boundary conditions 
%% default discretization is 0....M, with M a parameter
%% (in MATLAB the numbering is 1..M+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% example: 
%% fd1d_heat(5,0,1,0,0,1,0.01)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1, 
    M = 5; 
    a = 0; b = 1;
    ua = 0;
    ub = 0;
    T = 1;
    dt = .01;
end
clf;


dx = (b-a)/M;           %% use uniform grid
x = (a:dx:b)';

n = floor(T/dt);
t = dt:dt:T;

%%%%%%%%%%%%%%%%%%%%%%%% minus discrete laplacian
A = sparse(M+1,M+1);     
for i=2:M
    %% set up i'th row of the matrix
    A(i,i) = 2;
    A(i,i-1) = -1; 
    A(i,i+1) = -1; 
end
A =  A / (dx*dx);

uimpl = initfun(x);
Bimpl = sparse( eye(M+1,M+1) + dt*A );

%% apply/enforce Dirichlet boundary conditions
Bimpl (1,1) = 1;
Bimpl (M+1,M+1) = 1; 
fimpl (1,1) = ua;
fimpl (M+1,1) = ub;
    
maxerr = 0;
%% set up the matrix and the right hand side vector
for nt = 1:n
    t = nt*dt;
    uprev = uimpl;    
    
    fimpl (2:M,1) = dt*rhsfun(x(2:M),t) + uprev(2:M); 
    
    %% solve the linear system;
    uimpl = Bimpl \ fimpl;
    
    u = uimpl;
    exvec (1:M+1,1) = exfun(x(1:M+1),t);

    %% compute the error vector
    error = abs(exvec - u);
    %% report the error in max norm and in two discrete norm 
    maxerr = max(maxerr,norm(error,inf));
    
    %% for efficiency, you can comment out plotting, error report and pause   
    plot(x,exvec,'k',x,uimpl,'r*');axis([0 1 0 1]);
    fprintf('\nt = %g %g ',t,norm(error,inf)); 
    if nt> 4 pause(0.05);else pause; end;    
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = rhsfun(x,t)
v = zeros(size(x));

function v = initfun(x)
v = sin(pi*x);

function v = exfun(x,t)
v = sin(pi*x)*exp(-pi*pi*t);


    
