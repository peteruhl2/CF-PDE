%%% trying to get a traveling wave solution in a simpler CF PDE model. The
%%% model is
%%%
%%% f' = beta(1 - w/(1+w)) f (1-f) - d2 f - qfw + Df Lap(f)
%%% w' = lambda - mu w + Dw Lap(w)
%%%
%%% started 10/3/22

global r d lambda mu q c Df Dw

r = 5.0;
d = 0.1;
lambda = 1.0;
mu = 0.9;
q = 0.1;

c = 1.0;
Df = 1e-2;
Dw = 1e-0;

L = 5.0;
y0 = [0.5; 0.0; 1; 0.0];
domain = [-L L];

[t,y] = ode15s(@(t,y) rhs(t,y), domain, y0);

%%% Plots =================================================================

figure()
hold on; box on
plot(t,(y(:,1)),'linewidth',2);
plot(t,(y(:,3)),'linewidth',2);
legend('\phi_f','\phi_w')

figure()
hold on; box on
plot(t,y(:,2),'linewidth',2);
plot(t,y(:,4),'linewidth',2);
legend('\psi_f','\psi_w')






%%% Functions =============================================================

%%% ode function
function yp = rhs(t,y)
global r d lambda mu q c Df Dw

phif = y(1);
psif = y(2);
phiw = y(3);
psiw = y(4);

yp(1) = psif;
yp(2) = (-c*psif - r*(1 - phiw/(1+phiw))*phif*(1-phif) + d*phif + q*phif*phiw)/Df;
yp(3) = psiw;
% yp(4) = (-c*psiw - lambda + mu*phiw)/Dw;
yp(4) = (-c*psiw)/Dw;

yp = yp';
end