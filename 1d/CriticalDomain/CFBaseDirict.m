%%% CF PDE with Dirichlet BCs f is zero at the ends, oxygen is > 0 at the
%%% ends, not sure about c, maybe no flux
%%%
%%% started 9/22/2022

global beta dc df mu eta lambda Dc Df Dw L q

%%% assign variables
beta = 5.0;
dc = 1.0e-1;
df = 1.0e-1;

lambda = 0.0;
mu = 0.01;
eta = .50;
% eta = 2.14;
q = .5;

Dc = 4e-6;
Df = 4e-4;
Dw = 0.5e-1;
% Dw = 8.65;

%%% Domain
L = .800123;
x = linspace(-L,L,100);

% tolerance for finding radius
tol = 1e-2;

tmax = 100;
t = linspace(0,tmax);
dt = tmax/(length(t));

tic
m = 0;
sol = pdepe(m, @rhs, @icfun, @bcs, x, t);
toc

%%% Plots =================================================================

% % ---- Plot the solution --- %
% step=1;  %step every second
% for tt=1:step:length(t)
% % disp(sprintf('time = %d',(tt-1)*dt));
% disp((sprintf('time = %.02f out of %d',tt*dt,tmax)));
% plot(x,sol(tt,:,1),x,sol(tt,:,2),x,sol(tt,:,3),'linewidth',2);
% % legend('c','f','w');
% % ylim([-0.1,lambda/mu]);
% ylim([-0.1,1.1]);
% % axis tight
% drawnow;
% %     pause(0.0025);
% end
% close
% % -------------------------- %



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
global beta dc df mu eta Dc Df Dw L lambda q 

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2) - q*u(2)*u(3);
     - eta*u(1)*u(3)];

end

%%% IC function
function u0 = icfun(x)

% u0 = exp(-(x).^2);

u0 = [0.4*exp(-(x).^2); 
      0.2*exp(-(x).^2); 
      1 - exp(-(x).^2)];

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

