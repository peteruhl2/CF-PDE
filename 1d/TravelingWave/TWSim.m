%%% Implementing travling wave solution of CF PDE
%%% this program will solve the ode system to find the potential TW
%%% solution then use that as an IC in the PDE
%%%
%%% started 10/11/2022

global beta dc df mu eta local_lambda Dc Df Dw L q c
global c0 f0 w0 x

beta = 50.0;
dc = 1.0e-1;
df = 4.0e-0;

local_lambda = .10;
mu = 0.1;
eta = 1.0;
q = 5e-1;

%%% diffusion coefficients and wave speed
Dc = 1e-0;
Df = 1e-0;
Dw = 10e-4;
c = 10000;

%%% Domain
L = 50;
x = linspace(-L,L,100);

%%% solve ode for TW solution and interpolate solution
y0 = [.04; 0.; 2; 0.; local_lambda/mu; 0.];
domain = [-L L];
[t,y] = ode15s(@(t,y) TWrhs(t,y), domain, y0);

c0 = interp1(t,y(:,1),x);
f0 = interp1(t,y(:,3),x);
w0 = interp1(t,y(:,5),x);



%%% solve pde
tmax = 5;
t = linspace(0,tmax);
dt = tmax/(length(t));

tic
m = 0;
sol = pdepe(m, @rhs, @icfun, @bcs, x, t);
toc

%%% Plots =================================================================

% ---- Plot the solution --- %
step=1;  %step every second
for tt=1:step:length(t)
% disp(sprintf('time = %d',(tt-1)*dt));
disp((sprintf('time = %.02f out of %d',tt*dt,tmax)));
plot(x,sol(tt,:,1),x,sol(tt,:,2),x,sol(tt,:,3));
% legend('c','f','w');
% ylim([-0.1,lambda/mu]);
ylim([-0.1,1.1]);
% axis tight
drawnow;
%     pause(0.0025);
end
close
% -------------------------- %


[X,T] = meshgrid(x,t);

%%% aerobes
figure()
contourf(X,T,sol(:,:,1))
xlabel('x')
ylabel('t')
colorbar
title('Aerobes')

%%% anaerobes
figure()
contourf(X,T,sol(:,:,2))
xlabel('x')
ylabel('t')
colorbar
title('Anaerobes')

%%% oxygen
figure()
contourf(X,T,sol(:,:,3))
xlabel('x')
ylabel('t')
colorbar
title('Oxygen')











%%% Functions =============================================================

%%% TW ode function
function yp = TWrhs(t,y)
global beta dc df q local_lambda mu eta Dc Df Dw c

phic = y(1);
psic = y(2);
phif = y(3);
psif = y(4);
phiw = y(5);
psiw = y(6);

yp(1) = psic;
yp(2) = (-c*psic - (beta*phiw/(1+phiw))*phic*(1-phic-phif) + dc*phic)/Dc;
yp(3) = psif;
yp(4) = (-c*psif - beta*(1 - phiw/(1+phiw))*phif*(1-phic-phif) + df*phif + q*phif*phiw)/Df;
yp(5) = psiw;
yp(6) = (-c*psiw - local_lambda + mu*phiw + eta*phic*phiw)/Dw;

yp = yp';
end

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx)
global beta dc df mu eta Dc Df Dw L local_lambda q

% %%% assign spatial lambda value
% if abs(x) > 1.0*L
%     l = local_lambda;
% else
%     l = 0;
% end  

l = local_lambda;

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2) - q*u(2)*u(3);
     l - mu*u(3) - eta*u(1)*u(3)];

end

%%% IC function
function u0 = icfun(local_x)
global c0 f0 w0 x

% u0 = [0.4*exp(-(x).^2); 
%       0.2*exp(-(x).^2); 
%       1 - exp(-(x).^2)];

u0 = [interp1(x,c0,local_x); 
      interp1(x,f0,local_x);; 
      interp1(x,w0,local_x);];

end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)

% %%% Diriclet ul = ur = 0 for c,f
% %%% oxygen is no flux
% pl = [ul(1); ul(2); 0];
% ql = [0; 0; 1];
% pr = [ur(1); ur(2); 0];
% qr = [0; 0; 1];

%%% No flux
pl = [0; 0; 0];
ql = [1; 1; 1];
pr = [0; 0; 0];
qr = [1; 1; 1];

end