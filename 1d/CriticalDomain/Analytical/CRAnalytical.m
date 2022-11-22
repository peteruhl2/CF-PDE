%%% try to get a solver for fisher's equation
%%%
%%% u_t = D*u_xx + r*u*(1 - u)
%%%
%%% started 9/22/2022

global r D C L a q

r = 1.00;
C = 5/sqrt(6);
D = 5;
q = 2.1;
a = 10.5;

%%% Domain
L = 3.31096922;
x = linspace(-L,L);

tmax = 100;
t = linspace(0,tmax);
dt = tmax/(length(t));

m = 0;
sol = pdepe(m, @rhs, @icfun, @bcs, x, t);

% crit = pi*sqrt(D/(q-r));
crit = pi*sqrt(D/(r))/2;

%%% Plots =================================================================

% ---- Plot the solution --- %
step=1;  %step every second
for tt=1:step:length(t)
% disp(sprintf('time = %d',(tt-1)*dt));
disp((sprintf('time = %.02f out of %d',tt*dt,tmax)));
plot(x,sol(tt,:));
% legend('c','f','w');
% ylim([-0.1,lambda/mu]);
ylim([-0.1,1.1]);
drawnow;
%     pause(0.0025);
% waitforbuttonpress
end
% -------------------------- %



[X,T] = meshgrid(x,t);

figure()
contourf(X,T,sol)
xlabel('x')
ylabel('t')
colorbar



crit







%%% Functions =============================================================

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx)
global r D a q L

c = 1;
f = D*dudx;

% s = r*(1 - 1/(1 + exp(a*x)))*u*(1-u) - q*exp(-a*x)*u;

s = r*u*(1-u) - q*(exp(a*(x-L)) + exp(-a*(x+L)))*u;

% s = r*u*(1-u);
end

%%% IC function
function u0 = icfun(x)
global C L

%%% gaussian IC
u0 = exp(-(x).^2);

% % traveling wave solution
% u0 = (1 + C*exp((x+30)/sqrt(6*1.0)))^-2;

end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)
global D

%%% Diriclet ul = ur = 0
pl = ul;
ql = 0;
pr = ur;
qr = 0;

% %%% No flux
% pl = 0;
% ql = 1;
% pr = 0;
% qr = 1;

end

