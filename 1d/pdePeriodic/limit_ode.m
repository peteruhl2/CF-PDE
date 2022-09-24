%%% cf system with limit cycles
%%% 7/12/22

global beta beta2 b n dc df lambda mu eta k q

beta = 2;
beta2 = 2;
b = 1;
n = 1;
dc = 0.9;
df = 0.9;

lambda = 30;
mu = 0.001;
eta = 0.054;
q = 0;

k = 100*100;

C = 50;
F = 10;  
X = lambda/mu;

y0 = [C, F, X];
tspan = [0 5000];
[t, y] = ode15s(@(t,y) cf_eqs(t,y), tspan, y0);



figure()
hold on; box on
plot(t,y(:,1),'Linewidth',2)
plot(t,y(:,2),'Linewidth',2)
xlabel('Time')



%%% cf ode function
function yp = cf_eqs(t,y)
global beta beta2 b n dc df lambda mu eta k q 

c = y(1);
f = y(2);
x = y(3);

yp = zeros(3,1);

xx = x/k;

%    Cnext  = C + step*( rr1*C*(1-(C+F)/k)-d1*C )  
%    Fnext  = F + step*( rr2*F*(1-(C+F)/k) - d2*F -  q*F*X/k  )
%    Xnext  = X + step*( lambda - a*X- eta*X*C/k )

yp(1) = (beta*xx^n/(b^n + xx^n))*c*(1 - (c + f)/k) - dc*c;
yp(2) = (beta2*(1 - xx^n/(b^n + xx^n)))*f*(1 - (f + c)/k) - df*f - q*f*x/k;
yp(3) = lambda - mu*x - eta*(c)*x/k;


end