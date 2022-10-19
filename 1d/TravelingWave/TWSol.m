%%% first take at a traveling wave solution to the CF pde system
%%% model is
%%%
%%% c' = beta w/(1+w) c (1-c-f) - d1 c + Dc Lap(c)
%%% f' = beta(1 - w/(1+w)) f (1-c-f) - d2 f - qfw + Df Lap(f)
%%% w' = lambda - mu w - eta c w + Dw Lap(w)
%%%
%%% started 10/3/22

global beta dc df q lambda mu eta Dc Df Dw c

%%% params 
beta = 50.0;
b = 1.;
n = 1.0;
dc = 1.0e-1;
df = 4.0e0;
q = .5e1;
mu = .1;
lambda = .10;
eta = 10.0;

%%% diffusion coefficients and wave speed
Dc = 1e-0;
Df = 1e0;
Dw = 10e1;
c = 100;

L = .50;
y0 = [.04; 0.; 2; 0.; lambda/mu; 0.];
domain = [-L L];

[t,y] = ode15s(@(t,y) rhs(t,y), domain, y0);


%%% Plots =================================================================

figure()
hold on; box on
plot(t,(y(:,1)),'linewidth',2);
plot(t,(y(:,3)),'linewidth',2);
plot(t,(y(:,5)),'linewidth',2);
legend('\phi_c','\phi_f','\phi_w')

figure()
hold on; box on
plot(t,y(:,2),'linewidth',2);
plot(t,y(:,4),'linewidth',2);
plot(t,y(:,6),'linewidth',2);
legend('\psi_c','\psi_f','\psi_w')




%%% Functions =============================================================

%%% ode function
function yp = rhs(t,y)
global beta dc df q lambda mu eta Dc Df Dw c

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
yp(6) = (-c*psiw - lambda + mu*phiw + eta*phic*phiw)/Dw;

yp = yp';
end