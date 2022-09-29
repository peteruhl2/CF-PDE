%%% Critical domain size program. Checks if F survives for a range of
%%% domain sizes
%%% this one records the radius of the f blob after a long time (steady
%%% state-ish)
%%%
%%% started 9/27/22

% global beta dc df mu eta local_lambda Dc Df Dw L q

%%% assign variables
beta = 5.0;
dc = 1.0e-1;
df = 8.0e-1;

local_lambda = 10.0;
mu = 0.01;
eta = 5.0;
q = 4.5;

Dc = 4e-6;
Df = 4e-4;
Dw = 1.0e-3;

% last spot is for
p = [beta, dc, df, local_lambda, mu, eta, q, Dc, Df, Dw];

tmax = 30;
t = linspace(0,tmax,50);
dt = tmax/(length(t));

%%% Solver loop
DomainSize = linspace(0.1, 20, 100);
results = zeros(length(DomainSize),1);
m = 0;

% tolerance for finding radius
tol = 1e-2;

tic
parfor i = 1:length(results)
    i
    
    %%% Set domain
    L = DomainSize(i);
    x = linspace(-L,L,100);
    
    %%% solve pde
%     sol = pdepe(m, @rhs, @icfun, @bcs, x, t);

    fun = @(x,t,u,dudx) rhs(x,t,u,dudx,p,L);
    sol = pdepe(m, fun, @icfun, @bcs, x, t);
    
    %%% record radius of f blob
    thing = sol(end,:,2);
    idxs = find(thing > tol & x > 0);
    idx = max(idxs);
    
    %%% store 0 if there is no blob
    try
        results(i) = x(idx);
    catch
        results(i) = 0;
    end
    
end
toc

%%% Plots =================================================================

% figure()
hold on; box on
plot(DomainSize,smooth(results),'linewidth',2)
xlabel('Radius of Domain')
ylabel('Raidus of Anaerobe Blob')












%%% Functions =============================================================

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx,p,L)
% global beta dc df mu eta Dc Df Dw L local_lambda q

beta = p(1);
dc = p(2);
df = p(3);
local_lambda = p(4);
mu = p(5);
eta = p(6);
q = p(7);
Dc = p(8);
Df = p(9);
Dw = p(10);

%%% assign spatial lambda value
if abs(x) > 0.9*L
    l = local_lambda;
else
    l = 0;
end  

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2) - q*u(2)*u(3);
     l - mu*u(3) - eta*u(1)*u(3)];

end

%%% IC function
function u0 = icfun(x)

% u0 = exp(-(x).^2);

u0 = [0.4*exp(-(x).^2); 
      0.2*exp(-(x).^2); 
      1 - exp(-(x).^2)];


end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)

% %%% Diriclet ul = ur = 0 for c,f
% %%% oxygen is no flux
% pl = [ul(1); ul(2); 0];
% ql = [0; 0; 1];
% pr = [ur(1); ur(2); 0];
% qr = [0; 0; 1];

%%% No flux
pl = [0; 0; 0];
ql = [1; 1; 1];
pr = [0; 0; 0];
qr = [1; 1; 1];

end

%%% Spatial lambda function, return constant outside of a radius
function l = Lambda(x,L)
global local_lambda L

if abs(x) > 0.9*L
    l = local_lambda;
else
    l = 0;
end    

end