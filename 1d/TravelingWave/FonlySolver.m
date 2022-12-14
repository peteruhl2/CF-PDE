%%% F and oxygen only PDE model, just a solver
%%% model is
%%%
%%% f' = beta(1 - w/(1+w)) f (1-f) - d2 f - qfw + Df Lap(f)
%%% w' = lambda - mu w + Dw Lap(w)
%%%
%%% oxygen will have directlet at -L and L, > 0 on the right, 0 on the left
%%% F we'll see
%%%
%%% started 10/3/22

global r d lambda mu q c Df Dw
global x f0 w0 L

r = 15.0;
d = 1.5;
lambda = 3.0;
mu = 0.9;
q = 3;

c = 100000.0;
Df = 1e-4;
Dw = 1e-0;

%%% Domain
L = 3.;
x = linspace(-L,L,100);
% x = linspace(0, L ,100);

%%% solve ode for traveling wave solution
y0 = [0.1; 0.; lambda/mu; 0.];
domain = [-L L];

[t,y] = ode15s(@(t,y) TWrhs(t,y), domain, y0);

f0 = interp1(t,y(:,1),x);
w0 = interp1(t,y(:,3),x);

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
plot(x,sol(tt,:,1),x,sol(tt,:,2),'linewidth',2);
% legend('c','f','w');
% ylim([-0.1,lambda/mu]);
% ylim([-0.1,1.1]);
% axis tight
legend('Anaerobes','Oxygen')
drawnow;
%     pause(0.0025);
end
close
% -------------------------- %



[X,T] = meshgrid(x,t);

%%% anaerobes
figure()
contourf(X,T,sol(:,:,1))
xlabel('x')
ylabel('t')
colorbar
title('Anaerobes')

%%% oxygen
figure()
contourf(X,T,sol(:,:,2))
xlabel('x')
ylabel('t')
colorbar
title('Oxygen')


%%% Functions =============================================================

%%% ode function
function yp = TWrhs(t,y)
global r d lambda mu q c Df Dw

phif = y(1);
psif = y(2);
phiw = y(3);
psiw = y(4);

yp(1) = psif;
yp(2) = (-c*psif - r*(1 - phiw/(1+phiw))*phif*(1-phif) + d*phif + q*phif*phiw)/Df;
yp(3) = psiw;
yp(4) = (-c*psiw - lambda + mu*phiw)/Dw;

yp = yp';
end

%%% RHS PDE function
function [c_pdefunc,f,s] = rhs(x,t,u,dudx)
global r d lambda mu q Df Dw

c_pdefunc = [1; 1];
f = [Df*dudx(1); 
     Dw*dudx(2)];

s = [(r*(1 - u(2)/(1+u(2))))*u(1)*(1 - u(1)) - d*u(1) - q*u(1)*u(2);
     0*lambda - mu*u(2)];

end

%%% IC function
function u0 = icfun(local_x)

global x f0 w0 L

% u0 = [0.2*exp(-(x + L*0.7).^2); 
%       exp(-(x + L).^2)];

u0 = [interp1(x,f0,local_x);
      interp1(x,w0,local_x)];


% %%% exponential distriubtions at
% %%% left for ox and right for f
% u0 = [exp(-25*(local_x));
%       exp(25*(local_x - L))];
end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)
global lambda mu

%%% Diriclet ul = ur = 0
pl = [ul(1)*0; ul(2) - lambda/mu];
ql = [1; 0];
pr = [ur(1)*0; ur(2)];
qr = [1; 0];

% %%% No flux
% pl = [0; 0; 0];
% ql = [1; 1; 1];
% pr = [0; 0; 0];
% qr = [1; 1; 1];

end