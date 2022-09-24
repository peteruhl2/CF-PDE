%%% CF PDE with non-periodic boundaries
%%%
%%% started 9/22/2022

global r D C local_lambda

r = 10.00;
C = 5/sqrt(6);
D = 1.0;

local_lambda = 1.0;

%%% Domain
L = .7;
x = linspace(-L,L);

tmax = 20;
t = linspace(0,tmax);
dt = tmax/(length(t));

m = 0;
sol = pdepe(m, @rhs, @icfun, @bcs, x, t);

blah = zeros(length(x),1);
for i = 1:length(x)
    blah(i) = Lambda(x(i),L);
end


return

%%% Plots =================================================================

% % ---- Plot the solution --- %
% step=1;  %step every second
% for tt=1:step:length(t)
% % disp(sprintf('time = %d',(tt-1)*dt));
% disp((sprintf('time = %.02f out of %d',tt*dt,tmax)));
% plot(x,sol(tt,:));
% % legend('c','f','w');
% % ylim([-0.1,lambda/mu]);
% ylim([-0.1,1.1]);
% drawnow;
% %     pause(0.0025);
% end
% % -------------------------- %



[X,T] = meshgrid(x,t);

figure()
contourf(X,T,sol)
xlabel('x')
ylabel('t')
colorbar











%%% Functions =============================================================

%%% RHS PDE function
function [c,f,s] = rhs(x,t,u,dudx)
global r D

c = 1;
f = D*dudx;
s = r*u*(1-u);
end

%%% IC function
function u0 = icfun(x)
global C


u0 = exp(-(x).^2);

% % traveling wave solution
% u0 = (1 + C*exp((x+80)/sqrt(6*1.0)))^-2;

end

%%% BC function
function [pl,ql,pr,qr] = bcs(xl,ul,xr,ur,t)
global D

%%% Diriclet ul = ur = 0
pl = ul;
ql = 0;
pr = ur;
qr = 0;

% %%% No flux
% pl = 0;
% ql = 1;
% pr = 0;
% qr = 1;

end

%%% Spatial lambda function, return constant outside of a radius
function l = Lambda(x,L)
global local_lambda

if abs(x) > 0.95*L
    l = local_lambda;
else
    l = 0;
end    

end