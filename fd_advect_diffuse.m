function [btc,tsteps,nsol] = fd_advect_diffuse (filename,dt,ifshow,t1,t2,restart)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%% fd_advect_diffuse('m4test.fd',1e2,1,5)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% must be preceded by 
%% global mycase;mycase=1;global dirichlet_value;dirichlet_value=1;
%% genfiss(2,2,'m5fiss_300')
%% fd('m5fiss_300.fd');
%% generate_diffusion_dispersion_data
%% fd_advect_diffuse('m5fiss_300.fd',1,1,0,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% or for fluent data must be preceded by fluent data 
%%fluent2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DATA needed:
%% grid info: file of arbitrary name
%% diff11.dat, diff22.dat
%% adv_diff_bdary
%% vx.dat, vy,dat  and uvx.dat uvy.dat
%% porosity.dat
%% NNED globals: mycase, dirichlet_value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    t1 = 0; t2 = 10000;restart = 0;
end
if nargin < 6
    restart = 0;
end
if nargin < 3
    ifshow = 0;
end

%% control parameters
advect  = 1;
diffuse = 0;

if nargin < 3
    flag = 1;
end;

fid = fopen(filename,'r');
nx = fscanf(fid,'%d',1);
ny = fscanf(fid,'%d',1);
nz = fscanf(fid,'%d',1);

dx = fscanf(fid,'%g',nx);
dy = fscanf(fid,'%g',ny);
dz = fscanf(fid,'%g',nz);

x0 = fscanf(fid,'%g',1);
y0 = fscanf(fid,'%g',1);
z0 = fscanf(fid,'%g',1);
fclose(fid);

%%%%%%% diffusivivites
if diffuse > 0 
    permx=load('diff11.dat');
    permy=load('diff22.dat');
    permz=zeros(size(permx));

    %% load boundary information
    bfile = 'adv_diff_bdary.dat';
    fid = fopen(bfile,'r');
    ndirichlet = fscanf(fid,'%d',1);
    dirnodes=[];
    for i=1:ndirichlet  
        nod2 = fscanf(fid,'%g',2);
        face = fscanf(fid,'%g',1);
        dirnodes = [ dirnodes; [nod2(:,1)', face]];
    end   
    nneumann = fscanf(fid,'%d',1);
    neumnodes=[];
    for i=1:nneumann  
        nod2 = fscanf(fid,'%g',2);
        face = fscanf(fid,'%g',1);
        neumnodes = [ neumnodes; [nod2(:,1)', face]];
    end   
    fclose(fid);
end

%%%%%%%% load velocities
vx = load('vx.dat');[dum,ny] = size(vx);
vy = load('vy.dat');[nx,dum] = size(vy);
nz = 1;
%vz = load('vz.dat');
vz = zeros(nx,ny,2);

%max(vx(1,:,1)),pause
%max(vx(nx+1,:,1)),pause


%% check CFL condition
uvx = load('uvx.dat');
uvy = load('uvy.dat');
umag = sqrt(uvx.^2 + uvy.^2);
umax=max(umag(1:end));

if diffuse > 0 && advect > 0 
    pecletx=abs(umag./permx);
    peclety=abs(umag./permy);
    fprintf('Peclet number=%g umax=%g\n',max(max(pecletx(1:end)),max(peclety(1:end))),umax);
end


%% porosity and choice of time step honoring CFl condition 
por = load('porosity.dat');

%%%% input porosity - convert to pore volume
phi = por; %% pore volume 
for i=1:nx for j=1:ny for k=1:nz
        phi(i,j,k) = phi(i,j,k)*dx(i)*dy(j)*dz(k);  
    end;end;end;    

if advect > 0
    datcfl = phi./umag;
    dtcfl = min(datcfl(1:end));
    %dtcfl,max(umag(1:end)),pause
    
    dt0 = dtcfl;
    fprintf('dt chosen=%g  dtrequested=%g\n',dt0,dt);
    dt= min(dt0,dt);
end

tsteps=t1:dt:t2;
nt = length(tsteps);

%%%% create data structures, matrix and rhs 
x = dx; y = dy; z = dz ;
x(1)=x0; for i=2:nx x(i)=x(i-1)+dx(i);end;  xc = x + dx/2;
y(1)=y0; for i=2:ny y(i)=y(i-1)+dy(i);end;  yc = y + dy/2;
z(1)=z0; for i=2:nz z(i)=z(i-1)+dz(i);end;  zc = z + dz/2;

[xp,yp] = meshgrid(xc,yc);  xp = xp';  yp = yp';

if diffuse > 0 
[tx,ty,tz] = trans(dx,dy,dz,permx,permy,permz,dirnodes);

periodic = 0;
[nxyz,stiff,q] = formmatrhs(xc,yc,zc,dx,dy,dz,tx,ty,tz,...
        dirnodes,neumnodes,periodic,0);
stiff = dt.*stiff;
for m=1:nx for j=1:ny for k=1:nz
                nn = index(m,j,k,nx,ny,nz);
                stiff(nn,nn) = stiff(nn,nn) + phi(m,j,k);
    end;end;end;    
end
%%%%%% initial condition

if restart ~= 0
    nsolprev = load('concentration.dat');
    fprintf('Restarted !\n');
else
    for i=1:nx for j=1:ny for k=1:nz
        nsolprev(i,j,k) = 0;%1;%initfun ( xc(i), yc(j), zc(k));   
    end;end;end
end
clf;
btc = zeros(nt,1); 
for n = 1:nt
    
    t = t1 + n*dt;
    tsteps (n) = t;
    tprev = t-dt;
   
    total = nsolprev.*phi;
    totalc = sum(total(1:end));
   
    
    if advect
        concnew = advect_step(tprev,t,vx,vy,vz,phi,nsolprev);
        nsolprev = concnew;
    end
    
    if diffuse
%        max(nsolprev(1:end)),pause
        [nxyz,q] = formrhs(xc,yc,zc,dx,dy,dz,tx,ty,tz,...
            dirnodes,neumnodes,t);
        q = dt*q;
        for m=1:nx for j=1:ny for k=1:nz
                nn = index(m,j,k,nx,ny,nz);
                q (nn) = q (nn) + phi(m,j,k)*nsolprev(m,j,k);    
        end;end;end;    

        %%%%% solve the linear system
        sollin = stiff \ q;
        nsol = zeros(nx,ny,nz);
        %%%%%% unpack the solution
        for i=1:nx for j=1:ny for k=1:nz 
                nsol(i,j,k) = sollin(index(i,j,k,nx,ny,nz));    
            end; end; end;
    else
            nsol = nsolprev;
    end

    %% post-process and display solution 
        %fprintf('t=%g\n',t);
    btc(n) = sum(nsol(nx,:,1))./ny;    

    if ifshow > 0 && (mod(n,ifshow) == 0 || n==1 || n==nt)
            showtran(xc,yc,zc,xp,yp,nsol(:,:,1),btc,tsteps,n);pause(0.05);
             fprintf('t=%g nt=%d totalc=%g\n',t,n,totalc);
    end
    if ifshow < 0 && (mod(n,abs(ifshow)) == 0 || n==1 || n==nt)
     %       showtran(xc,yc,zc,xp,yp,nsol(:,:,1),btc,tsteps,n);pause(0.05);
             fprintf('t=%g nt=%d totalc=%g\n',t,n,totalc);
    end
    nsolprev = nsol;
end   %% time stepping loop 

save('concentration.dat','nsol','-ASCII','-DOUBLE','-TABS');
save('xc.dat','xc','-ASCII','-DOUBLE','-TABS');
save('yc.dat','yc','-ASCII','-DOUBLE','-TABS');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% actual advection time step
function conc =  advect_step(told,tnew,vx,vy,vz,phi,concold)
inflow = 0;
[nx,ny,nz] = size(concold);
conc = concold.*phi;
dt = tnew - told;
t = tnew;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% advect in x - direction
for j=1:ny 
    for k=1:nz 
        for l=1:nx+1 
    quant=0;
%%  find the quantity to be advected, should take into account boundary
%%  conditions
   if l>=2 && l<= nx
     if vx(l,j,k) >= 0
       quant = concold(l-1,j,k);
     else
       quant = concold(l,j,k);
     end;
   elseif l == 1
     if vx(l,j,k) >= 0
        dirval = dirichlet_fun (xc(l) - dx(l)/2, yc(j), zc(k),t) ; 
        quant = dirval;
%        quant,pause
     else
       quant = concold(l,j,k);
     end;
   else % i == nx+1
     if vx(l,j,k) >= 0
       quant = concold(l-1,j,k);
     else
         dirval = dirichlet_fun (xc(l-1) + dx(l-1)/2, yc(j), zc(k),t) ; 
       quant = dirval;
     end;
   end;
  
   if l<=nx 
       conc(l,j,k)   = conc(l,j,k) + dt*vx (l,j,k)*quant;
      inflow = inflow + dt*vx(l,j,k)*quant;
      
   end;  %%% or outflux 
   if l>1  %% or boundary condition
       conc(l-1,j,k) = conc(l-1,j,k) -  dt*vx (l,j,k)*quant;
       inflow = inflow - dt*vx(l,j,k)*quant;
   end;
        end
    end
end
%inflow,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% advect in y - direction
for i=1:nx 
    for k=1:nz 
        for j=1:ny+1 

%%  find the quantity to be advected, should take into account boundary
%%  conditions
   if j>=2 && j<= ny
     if vy(i,j,k) >= 0
       quant = concold(i,j-1,k);
     else
       quant = concold(i,j,k);
     end;
   elseif j == 1
     if vy(i,j,k) >= 0
         dirval = dirichlet_fun (xc(i), yc(j) -dy(j)/2 , zc(k),t) ; 
       quant = dirval;
     else
       quant = concold(i,j,k);
     end;
   else % j == ny+1
     if vy(i,j,k) >= 0
       quant = concold(i,j-1,k);
     else
         dirval = dirichlet_fun (xc(i), yc(j-1) +dy(j-1)/2 , zc(k),t) ; 
       quant = dirval;
     end;
   end;
   
   if j<=ny
       conc(i,j,k)   = conc(i,j,k) + dt*vy (i,j,k)*quant;
       inflow = inflow + dt*vy (i,j,k)*quant;
   end;  %%% or outflux 
   if j>1  %% or boundary condition
       conc(i,j-1,k) = conc(i,j-1,k) -  dt*vy (i,j,k)*quant;
       inflow = inflow - dt*vy (i,j,k)*quant;
   end;
end;end;end;
%inflow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% advect in z - direction
for j=1:ny for i=1:nx for k=1:nz+1 

%%  find the quantity to be advected, should take into account boundary
%%  conditions
   if k>=2 && k<= nz
     if vz(i,j,k) >= 0
       quant = concold(i,j,k-1);
     else
       quant = concold(i,j,k);
     end;
   elseif k == 1
     if vz(i,j,k) >= 0
         dirval = dirichlet_fun (xc(i), yc(j), zc(k)-dz(k)/2,t) ; 
       quant = dirval;
     else
       quant = concold(i,j,k);
     end;
   else % k == nz+1
     if vz(i,j,k) >= 0
       quant = concold(i,j,k-1);
     else
         dirval = dirichlet_fun (xc(i), yc(j), zc(k-1)+dz(k-1)/2,t) ; 
       quant = dirval;
     end;
   end;
   
   if k<=nz 
       conc(i,j,k)   = conc(i,j,k) + dt*vz (i,j,k)*quant;
       inflow = inflow + dt*vz (i,j,k)*quant;
   end;  %%% or outflux 
   if k>1  %% or boundary condition
       conc(i,j,k-1) = conc(i,j,k-1) -  dt*vz (i,j,k)*quant;
       inflow = inflow - dt*vz (i,j,k)*quant;
   end;
end;end;end;
%inflow,

%%%%%%%%%% post-process: bring back the old values 
conc = conc ./phi;
%max(conc(1:end)),pause

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function post = postvec(vec,par)
[nx,ny,nz] = size(vec);

switch par
    case 1
        post = zeros(nx-1,ny,nz);
        for i=1:nx-1 post(i,:,:) = 1/2*(vec(i,:,:)+vec(i+1,:,:));end
    case 2
        post = zeros(nx,ny-1,nz);
        for i=1:ny-1 post(:,i,:) = 1/2*(vec(:,i,:)+vec(:,i+1,:));end
end

end

%%% contours of numerical and analytical solutions
function showtran(xc,yc,zc,xp,yp,conc,btc,tsteps,nt)

ifregular =0;

if ifregular
cross1=6;
cross2=12;

mytit=sprintf('Solution at t=%g',tsteps(nt)); 
mytit1=sprintf('Solution at t=%g y=%g',tsteps(nt),yc(cross1)); 
mytit2=sprintf('Solution at t=%g y=%g',tsteps(nt),yc(cross2)); 

subplot(1,4,1);plot(xc,conc(:,cross1,1));axis([0 1 0 1.2]);
title(mytit1);
subplot(1,4,2);plot(xc,conc(:,cross2,1));axis([0 1 0 1.2]);
title(mytit2);

subplot(1,4,3);
contourf(xp,yp,conc);axis square; caxis([0 1]);axis square; colorbar;

title(mytit);
%pcolor(conc(:,:,1)');colorbar;

%fprintf('%s max conc = %g\n',mytit,max(conc(1:end)));

subplot(1,4,4);semilogx(tsteps,btc);axis([0 tsteps(end) 0 1.2]);
title('Breakthrough curve');
else
%    figure(1);
    pcolor(conc');colorbar;
%    contourf(xp,yp,conc);axis square; caxis([0 1]);colorbar;
end


end