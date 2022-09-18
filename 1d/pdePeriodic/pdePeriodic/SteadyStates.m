%%% steady states for the model
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


%%% 1 extinction
c1 = 0;
f1 = 0;
w1 = lambda/mu;

%%% 2 c only
c2top = beta*lambda - d1*lambda - d1*mu;
c2bottom = d1*eta + beta*lambda;
c2 = c2top/c2bottom;

f2 = 0;

w2top = (-d1*eta - beta*lambda);
w2bottom = d1*eta - beta*eta - beta*mu;
w2 = w2top/w2bottom;

%%% 3 f only
c3 = 0;
f3 = (beta*mu - d2*lambda - d2*mu)/(beta*mu);
w3 = lambda/mu;

%%% 4 coexistence
c4 = (d2*lambda - d1*mu)/(d1*eta);
f4 = (d1*beta*eta - eta*d1^2 - d1*d2*eta - d2*beta*lambda + d1*beta*mu)/(d1*beta*eta);
w4 = d1/d2;



E1 = jac(c1,f1,w1);
E2 = jac(c2,f2,w2);
E3 = jac(c3,f3,w3);
E4 = jac(c4,f4,w4);

ss1 = eig(E1)
ss2 = eig(E2)
ss3 = eig(E3)
ss4 = eig(E4)









%%% Functions =============================================================

function J = jac(c,f,w)
global beta d1 d2 mu eta

J(1,1) = -d1 - (c*w*beta)/(1+w) + (w*beta*(1-c-f))/(1+w);
J(1,2) = -(c*w*beta)/(1+w);
J(1,3) = -(c*(1-c-f)*w*beta)/((1+w)^2) + (c*(1-c-f)*beta)/(1+w);

J(2,1) = -f*(1 - w/(1+w))*beta;
J(2,2) = -d2 + beta*(1-c-f)*(1-w/(1+w)) - beta*f*(1 - w/(1+w));
J(2,3) = beta*(1-c-f)*f*(w/((1+w)^2) - 1/(1+w));

J(3,1) = -eta*w;
J(3,2) = 0;
J(3,3) = -eta*c - mu;

end
