%%% CF PDE with Dirichlet BCs f is zero at the ends, oxygen is > 0 at the
%%% ends, not sure about c, maybe no flux
%%%
%%% started 9/22/2022

global beta dc df eta Dc Df Dw L q a

%%% assign variables
beta = 1.0;
dc = 1.0e-4;
df = 1.0e-4;

eta = 2.50;
q = 3.5;

Dc = 4e-8;
Df = 4e-8;
Dw = 1;

%%% oxygen decent parameter
a = 20;

%%% Domain
L = 5.000123;
x = linspace(-L,L,100);

% tolerance for finding radius
tol = 1e-2;

tmax = 400;
t = linspace(0,tmax,100);
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
plot(x,sol(tt,:,1),x,sol(tt,:,2),x,sol(tt,:,3),'linewidth',2);
% legend('c','f','w');
% ylim([-0.1,lambda/mu]);
ylim([-0.1,1.1]);
% axis tight
drawnow;
%     pause(0.0025);
% waitforbuttonpress
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

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx)
global beta dc df eta Dc Df Dw L q 

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2) - q*u(2)*u(3);
     - eta*u(1)*u(3)];

end

%%% IC function
function u0 = icfun(x)
global L a

u0 = [0.4*exp(-(x).^2); 
      0.4*exp(-(x).^2); 
      (exp(a*(x-L)) + exp(-a*(x+L)))];

end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)

%%% Diriclet ul = ur = 0 for oxygen
%%% c,f are no flux
pl = [0; 0; ul(3) - 1];
ql = [1; 1; 0];
pr = [0; 0; ur(3) - 1];
qr = [1; 1; 0];

% %%% No flux
% pl = [0; 0; 0];
% ql = [1; 1; 1];
% pr = [0; 0; 0];
% qr = [1; 1; 1];

end

