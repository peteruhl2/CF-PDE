%%% CF PDE with non-periodic boundaries
%%%
%%% started 9/22/2022

global beta dc df mu eta local_lambda Dc Df Dw L q

%%% assign variables
beta = 5.0;
dc = 1.0e-1;
df = 8.0e-1;

local_lambda = 10.0;
mu = 0.01;
eta = 5.0;
q = 4.5;

Dc = 4e-6;
Df = 4e-4;
Dw = 0.5e0;

%%% Domain
L = 0.97;
x = linspace(-L,L,100);

tmax = 480;
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
plot(x,sol(tt,:,1),x,sol(tt,:,2),x,sol(tt,:,3),'linewidth',2);
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

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx)
global beta dc df mu eta Dc Df Dw L local_lambda q

%%% assign spatial lambda value
if abs(x) > 0.9*L
    l = local_lambda;
else
    l = 0;
end  

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2) - q*u(2)*u(3);
     l - mu*u(3) - eta*u(1)*u(3)];

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

%%% Spatial lambda function, return constant outside of a radius
function l = Lambda(x,L)
global local_lambda L

if abs(x) > 0.9*L
    l = local_lambda;
else
    l = 0;
end    

end