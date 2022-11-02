%%% Critical domain size program. Gets the critical domain size for a range
%%% of eta (consumption) values
%%%
%%% started 11/1/22

% global beta dc df mu eta local_lambda Dc Df Dw L q

%%% assign variables
beta = 5.0;
dc = 1.0e-1;
df = 1.0e-1;

lambda = 0.0;
mu = 0.01;
eta = .50;
q = .5;

Dc = 4e-6;
Df = 4e-4;
Dw = 0.5e-1;

p = [beta, dc, df, q, Dc, Df, Dw];

tmax = 1080;
t = linspace(0,tmax,50);
dt = tmax/(length(t));

%%% Solver loop
res = 200;
Eta = linspace(0,0.5,res);
DomainSize = linspace(0.01, 5, res);
results = zeros(length(Eta),1);
tempresult = 0;
m = 0;

% tolerance for finding radius
tol = 1e-2;

tic
parfor i = 1:length(Eta)
    
    eta = Eta(i);
    tempresult = 0;
    
    for j = 1:length(DomainSize)    
        [i j]

        %%% Set domain
        L = DomainSize(j);
        x = linspace(-L,L,100);

        %%% solve pde
        fun = @(x,t,u,dudx) rhs(x,t,u,dudx,p,eta);
        sol = pdepe(m, fun, @icfun, @bcs, x, t);

        %%% record radius of f blob
        thing = sol(end,:,2);
        idxs = find(thing > tol & x > 0);
        idx = max(idxs);

        %%% store 0 if there is no blob
        try
            tempresult = x(idx);
        catch
            tempresult = 0;
        end
        
        %%% break if passed critical domain size
        if tempresult > 0
            %%% get the previous x point, which is 
            %%% the last one where F went extinct
            results(i) = L;
            break
        end
        
    end % while loop
end % loop on parameter loop
toc

% %%% print critical domain size
% critdx = find(results <= 0, 1, 'last' );
% critd = DomainSize(critdx)




%%% Plots =================================================================

figure()
hold on; box on
plot(Eta,(results),'linewidth',2)
xlabel('\eta')
ylabel('Critical Domain Size')












%%% Functions =============================================================
p = [beta, dc, df, q, Dc, Df, Dw];
%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx,p,eta)
% global beta dc df mu eta Dc Df Dw L local_lambda q

beta = p(1);
dc = p(2);
df = p(3);
q = p(4);
Dc = p(5);
Df = p(6);
Dw = p(7);

c = [1; 1; 1];
f = [Dc*dudx(1); Df*dudx(2); Dw*dudx(3)];

s = [(beta*u(3)/(1 + u(3)))*u(1)*(1 - u(1) - u(2)) - dc*u(1);
     (beta*(1 - u(3)/(1+u(3))))*u(2)*(1 - u(1) - u(2)) - df*u(2) - q*u(2)*u(3);
     - eta*u(1)*u(3)];

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

%%% Diriclet ul = ur = 0 for oxygen
%%% c,f are no flux
pl = [0; 0; ul(3) - 1];
ql = [1; 1; 0];
pr = [0; 0; ur(3) - 1];
qr = [1; 1; 0];

% %%% No flux
% pl = [0; 0; 0];
% ql = [1; 1; 1];
% pr = [0; 0; 0];
% qr = [1; 1; 1];

end