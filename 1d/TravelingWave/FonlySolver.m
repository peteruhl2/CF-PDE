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

r = 5.0;
d = 1.0;
lambda = 1.0;
mu = 0.9;
q = 0;

c = 10.0;
Df = 1e-2;
Dw = 1e-0;

%%% Domain
L = 20.;
x = linspace(-L,L,100);

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

%%% RHS PDE function
function [c_pdefunc,f,s] = rhs(x,t,u,dudx)
global r d lambda mu q c Df Dw

c_pdefunc = [1; 1];
f = [Df*dudx(1); Dw*dudx(2)];

s = [(r*(1 - u(2)/(1+u(2))))*u(1)*(1 - u(1)) - d*u(1) - q*u(1)*u(2);
     0*lambda - mu*u(2)];

end

%%% IC function
function u0 = icfun(x)

global L

u0 = [0.2*exp(-(x + L*0.7).^2); 
      exp(-(x + L).^2)];

end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)
global lambda mu

%%% Diriclet ul = ur = 0 for c,f
%%% oxygen is no flux
pl = [ul(1); ul(2) - lambda/mu];
ql = [0; 0];
pr = [ur(1); ur(2)];
qr = [0; 0];

% %%% No flux
% pl = [0; 0; 0];
% ql = [1; 1; 1];
% pr = [0; 0; 0];
% qr = [1; 1; 1];

end