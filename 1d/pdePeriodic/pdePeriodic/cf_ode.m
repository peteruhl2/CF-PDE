%%% simple CF ODE system to compare with PDE model
%%% started 8/29/22

close all

global beta b n dc df q mu lambda eta 

beta = .68;
b = 1.;
n = 1.0;
dc = 9e-3;
df = 9e-3;
q = 4.48e-0;
mu = 0.01312;
lambda = 0.01312;
eta = 2.54;

c0 = 0.4;
f0 = 0.3;
w0 = 0.1;
y0 = [c0, f0, w0];

tspan = [0, 500];

[t,y] = ode15s(@(t,y) rhs(t,y), tspan, y0);


% plots ===================================================================
figure()
hold on; box on
plot(t,y(:,1),'linewidth',2)
plot(t,y(:,2),'linewidth',2)
xlabel('Time')
legend('C','F','Fontsize',18,'Location','east')

figure()
hold on; box on;
plot(t, y(:,3), 'linewidth',2)
xlabel('Time')
legend('Oxygen','Fontsize',18)






% functions ===============================================================

function yp = rhs(t,y)
global beta b n dc df q mu lambda eta 

c = y(1);
f = y(2);
w = y(3);

yp(1) = (beta*w^n)/(b^n + w^n)*c*(1-c-f) - dc*c;
yp(2) = beta*(1 - (w^n)/(b^n + w^n))*f*(1-c-f) - df*f - q*f*w;
yp(3) = lambda - mu*w - eta*c*w;

% yp(1) = (beta*w^n)/(b^n + w^n)*c*(1-c-f) - dc*c;
% yp(2) = - df*f - q*f*w;
% yp(3) = lambda - mu*w - eta*c*w;


yp = yp';

end

