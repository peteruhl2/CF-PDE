%%% jacobian and steady-states for CF model with constant growth
%%% but toxicity and consumption terms
%%% started 9/9/22

global r1 r2 d1 d2 q mu lambda eta 

r1 = 2;
r2 = 2.4;
d1 = 9e-1;
d2 = 9e-1;
q = 4.48e-9;
mu = 0.051312;
lambda = 2.01312;
eta = .054;


%%% steady state solutions
ctop = d1*r2*mu - q*r1*lambda - d2*r1*mu;
cbottom = eta*(d2*r1 - d1*r2);
c = ctop/cbottom;

ftop = r2*eta*d1^2 + ( q*lambda+d2*( eta+mu ) )*r1^2 - d1*r1*( d2*eta+r2*( eta+mu ) );
fbottom = r1*eta*(d2*r1 - d1*r2);
f = ftop/fbottom;

w = (d1*r2 - d2*r1)/(q*r1);


%%% plot characteristic equation
AA = r1 - d2;
BB = r2 - d2;
xi = linspace(-5.6,5.2);

y = (-eta*w)*(r1*q*c*f) ... 
    - (eta*c - xi).*( (AA-2*r1*c - r1*f - xi) .*( BB-r2*c - 2*r2*f - q*w - xi) - r1*r2*c*f );


hold on
plot(xi,y,'linewidth',2)
yline(0)
xlabel('\xi')
ylabel('Char. Eq.')
ylim([-1,1])


% discrim(1,1,-3,1)

%%% Characteristic equation stuff =========================================

A = r1-2*r1*c - r1*f - d1;
B = r2 - r2*c - 2*r2*f - d2 - q*w;
E = r1*c*eta*w*q*f - r1*r2*eta*f*c^2;

bb = eta*c - A - B;
cc = A*B - eta*c*A - eta*c*B - r1*r2*c*f;
dd = E + eta*c*A*B;



delta = discrim(1,bb,cc,dd)

%%% Functions =============================================================

function d = discrim(a,b,c,d)

d = 18*a*b*c*d - 4*d*b^3 + (b*c)^2 - 4*a*c^3 - 27*(a*d)^2;

end











