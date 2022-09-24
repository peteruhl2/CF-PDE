%%% for periodic BCs use
% ipx=[Nx,1:Nx-1];
% imx=[2:Nx,1];
% ipy=[Ny,1:Ny-1];
% imy=[2:Ny,1];

% Physical parameters
sig=1;
myalpha=1;
tf=25;
Lx=10;
Ly=10;

%numerical constants
disp=100; % number of snapshots to display and save
Nx=100;
Ny=100;
dt=0.01;
x=linspace(-Lx,Lx,Nx);
y=linspace(-Ly,Ly,Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);
ootdx=1/(2*dx);
ootdy=1/(2*dy);
oodx2=1/(dx^2);
oody2=1/(dy^2);
ipx=[1,1:Nx-1];
imx=[2:Nx,Nx];
ipy=[1,1:Ny-1];
imy=[2:Ny,Ny];

% Define matrices to store solutions
allu=zeros(disp+1,Ny,Nx);
allt=zeros(1,disp+1);

% Initial condition
[X,Y]=meshgrid(x,y);
R2=X.^2+Y.^2;
u0=exp(-R2./(2*sig));
u=u0;

%plot initial condition
figure(1);clf
surf(x,y,u0)
xlabel('x')
ylabel('y')
zlabel('u(x,y,t=0)')
title('Initial condition')
axis tight;
ax=axis;

ii=0;
isave=0;
for t=0:dt:tf
 ii=ii+1;
 
 %compute uxx and uyy:   
 uxx = oodx2*(u(:,ipx)-2*u+u(:,imx)); 
 uyy = oody2*(u(ipy,:)-2*u+u(imy,:)); 
 
 %compute next time:
 uN = u + (myalpha*dt)*(uxx+uyy);
 
 %update:
 u=uN;
 
 %time to save the solution
 if(mod(ii-1,(tf/dt)/disp)==0)
  isave=isave+1;
  allt(isave)=t;
  allu(isave,:,:)=u;
  figure(2);
  surf(x,y,u)
  xlabel('x')
  ylabel('y')
  zlabel('u(x,y,t)')
  title(['t=',num2str(t),'(/',num2str(tf),')'])
  axis(ax)
  %pause(0.01)
  drawnow
 end
end



