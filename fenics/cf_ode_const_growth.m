%%% simple CF ODE system to compare with PDE model
%%% started 8/29/22

close all

global r1 r2 d1 d2 q mu lambda eta 

r1 = 2;
r2 = 2.2;
d1 = 9e-10;
d2 = 9e-1*0;
q = 4.48e-10;
mu = 0.01312;
lambda = .01312;
eta = .054;

c0 = 0.4;
f0 = 0.3;
w0 = 0.1;
y0 = [c0, f0, w0];

tspan = [0, 5000];

[t,y] = ode15s(@(t,y) rhs(t,y), tspan, y0);


% plots ===================================================================
figure()
hold on; box on
plot(t,y(:,1),'linewidth',2)
plot(t,y(:,2),'linewidth',2)
plot(t,y(:,3),'linewidth',2)
xlabel('Time')
legend('C','F','X','Fontsize',18,'Location','east')

% figure()
% hold on; box on;
% plot(t, y(:,3), 'linewidth',2)
% xlabel('Time')
% legend('Oxygen','Fontsize',18)



%%% steady state solutions
ctop = d1*r2*mu - q*r1*lambda - d2*r1*mu;
cbottom = eta*(d2*r1 - d1*r2);
cstar = ctop/cbottom

ftop = r2*eta*d1^2 + ( q*lambda+d2*( eta+mu ) )*r1^2 - d1*r1*( d2*eta+r2*( eta+mu ) );
fbottom = r1*eta*(d2*r1 - d1*r2);
fstar = ftop/fbottom

wstar = (d1*r2 - d2*r1)/(q*r1)


% functions ===============================================================

function yp = rhs(t,y)
global r1 r2 d1 d2 q mu lambda eta 

c = y(1);
f = y(2);
w = y(3);

yp(1) = r1*c*(1-c-f) - d1*c;
yp(2) = r2*f*(1-c-f) - d2*f - q*f*w;
yp(3) = lambda - mu*w - eta*c*w;

yp = yp';

end