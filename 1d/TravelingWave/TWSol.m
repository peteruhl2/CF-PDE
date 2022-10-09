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
beta = 5.0;
b = 1.;
n = 1.0;
dc = 1.0e-1;
df = 8.0e-1;
q = .5e-0;
mu = .1;
lambda = .10;
eta = 10.0;

%%% diffusion coefficients and wave speed
Dc = 1e-1;
Df = 1e-1;
Dw = 10e-0;
c = 10000;

y0 = [.4; 0.1; 0.3; 0.1; lambda/mu; 0.1];
tspan = [0 1];

[t,y] = ode15s(@(t,y) rhs(t,y), tspan, y0);


%%% Plots =================================================================

hold on; box on
plot(t,y(:,1),'linewidth',2);
plot(t,y(:,2),'linewidth',2);
plot(t,y(:,3),'linewidth',2);
plot(t,y(:,4),'linewidth',2);
plot(t,y(:,5),'linewidth',2);
plot(t,y(:,6),'linewidth',2);
legend('x_1','x_2','x_3','x_4','x_5','x_6')





%%% Functions =============================================================

%%% ode function
function yp = rhs(t,y)
global beta dc df q lambda mu eta Dc Df Dw c

x1 = y(1);
x2 = y(2);
x3 = y(3);
x4 = y(4);
x5 = y(5);
x6 = y(6);

yp(1) = (-c*x1 - (beta*x6/(1+x6))*x2*(1-x2-x4) + dc*x2)/Dc;
yp(2) = x1;
yp(3) = (-c*x3 - (beta*(1 - x6/(1+x6)))*x4*(1-x2-x4) + df*x4 + q*x4*x6)/Df;
yp(4) = x3;
yp(5) = (-c*x5 - lambda + mu*x6 + eta*x2*x6)/Dw;
yp(6) = x5;

yp = yp';
end