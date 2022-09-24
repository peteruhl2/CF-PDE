%%% solver for the model
%%% c' = beta*w/(1+w)*c*(1-c-f) - d1c
%%% f' = beta*(1 - w/(1+w))*f*(1-c-f) - d2f
%%% w' = lambda - mu*w - eta*c*w
%%%
%%% started 9/16/22

global beta d1 d2 lambda mu eta

beta = 2.0;
d1 = 0.1;
d2 = 0.015;
lambda = 2.6;
mu = 0.2;
eta = 1.1;

c0 = 0.4;
f0 = 0.3;
w0 = 1;

y0 = [c0, f0, w0];
tspan = [0 1000];
[t, y] = ode45(@(t,y) rhs(t,y), tspan, y0);

figure()
hold on
plot(t, y(:,1), 'Linewidth', 2)
plot(t, y(:,2), 'Linewidth', 2)
plot(t, y(:,3), 'Linewidth', 2)
legend('C','F','W')









%%% Functions =============================================================

function yp = rhs(t,y)
global beta d1 d2 lambda mu eta

c = y(1);
f = y(2);
w = y(3);

yp(1) = (beta*w/(1+w))*c*(1-c-f) - d1*c;
yp(2) = (beta*(1 - w/(1+w)))*f*(1-c-f) - d2*f;
yp(3) = lambda - mu*w - eta*c*w;

yp = yp';

end