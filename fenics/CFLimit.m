%%% simple CF ODE system with limit cycles
%%% started 9/10/22

% close all

global beta b n dc df q mu lambda eta k

k = 100*100;
beta = 2.0;
b = 1.;
n = 1.0;
dc = 9e-1;
df = 9e-1;
q = 4.48e-4;
mu = 0.001;
lambda = 0.002;
eta = 5.4e-6;


c0 = 100.4;
f0 = 2.3;
w0 = 1e-5;
y0 = [c0; f0; w0];

tspan = [0, 5000];

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

% hold on; box on
% plot(y(:,1),y(:,2),'linewidth',1)
% scatter(y(end,1),y(end,2))
% xlabel('C')
% ylabel('F')






% functions ===============================================================

function yp = rhs(t,y)
global beta b n dc df q mu lambda eta k

c = y(1);
f = y(2);
w = y(3);

yp(1) = (beta*w^n)/(b^n + w^n)*c*(1-(c+f)/k) - dc*c;
yp(2) = beta*(1 - (w^n)/(b^n + w^n))*f*(1-(c+f)/k) - df*f - q*f*w;
yp(3) = lambda - mu*w - eta*c*w;

% yp(1) = (beta*w^n)/(b^n + w^n)*c*(1-c-f) - dc*c;
% yp(2) = - df*f - q*f*w;
% yp(3) = lambda - mu*w - eta*c*w;


yp = yp';

end