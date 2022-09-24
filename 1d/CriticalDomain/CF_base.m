%%% CF PDE with non-periodic boundaries
%%%
%%% started 9/22/2022

global beta dc df mu eta local_lambda Dc Df Dw L

beta = 2.0;
dc = 0.01;
df = 0.01;

local_lambda = 1.0;
mu = 1;
eta = 0;

Dc = 4e-1;
Df = 4e-1;
Dw = 4e-1;

%%% Domain
L = 5.7;
x = linspace(-L,L);

tmax = 20;
t = linspace(0,tmax);
dt = tmax/(length(t));

m = 0;
sol = pdepe(m, @rhs, @icfun, @bcs, x, t);


%%% Plots =================================================================

% % ---- Plot the solution --- %
% step=1;  %step every second
% for tt=1:step:length(t)
% % disp(sprintf('time = %d',(tt-1)*dt));
% disp((sprintf('time = %.02f out of %d',tt*dt,tmax)));
% plot(x,sol(tt,:));
% % legend('c','f','w');
% % ylim([-0.1,lambda/mu]);
% ylim([-0.1,1.1]);
% drawnow;
% %     pause(0.0025);
% end
% % -------------------------- %



[X,T] = meshgrid(x,t);

figure()
contourf(X,T,sol(:,:,3))
xlabel('x')
ylabel('t')
colorbar
title('Oxygen')











%%% Functions =============================================================

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx)
global beta dc df mu eta Dc Df Dw L local_lambda

%%% assign spatial lambda value
if abs(x) > 0.9*L
    l = local_lambda;
else
    l = 0;
end  

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2);
     l - mu*u(3) - eta*u(1)*u(3)];
end

%%% IC function
function u0 = icfun(x)

% u0 = exp(-(x).^2);

u0 = [1; 1; exp(-(x).^2)];


end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)

%%% Diriclet ul = ur = 0 for c,f
%%% oxygen is no flux
pl = [ul(1); ul(2); 0];
ql = [0; 0; 1];
pr = [ur(1); ur(2); 0];
qr = [0; 0; 1];

% %%% No flux
% pl = 0;
% ql = 1;
% pr = 0;
% qr = 1;

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