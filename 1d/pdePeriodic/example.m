%%% Example from a paper

global ru rv rw bu bv su sw ss Du Dv Dw

ru = 1.0;
rv = 1.0;
rw = 0.01;
bu = 0.2;
bv = 0;
su = 0;
sw = 2.0;
ss = 1.0;


Du = 1.32e-0;
Dv = 10.32e-0;
Dw = 100.32e-0;

% Dc = 0; Df = 0; Dw = 0;

xlist = linspace(0,1,100);
nx = length(xlist);

dt = 1;
tlist=0:dt:5000;

% ic = [0.600*ones(1,nx);
%       0.6*ones(1,nx);
%       (0.6)*ones(1,nx)];
  
%%% maybe do gausians
% ic(1,55) = 0.2;
% ic(1,10) = 0.1;

ic(1,:) = ic(1,:) + 0.6*exp(-(100*(xlist - 0.4).^2)/0.1);
% ic(1,:) = ic(1,:) + 0.1*exp(-(10*(xlist - 0.7).^2)/0.01);
% ic(1,:) = ic(1,:) + 0.1*exp(-(10*(xlist - 0.9).^2)/0.01);

ic(2,:) = 0.5*exp(-(15*(xlist - 0.6).^2)/0.05);


% ic(2,48) = 0.3;
% ic(2,58) = 0.4;
% ic(2,10) = 0.4;


%%% Solve
disp('solving pde');
sol = pbcpdeSolver(@my_pde,ic,xlist,tlist);
disp('done')

% % ---- Plot the solution --- %
% step=1/dt;  %step every second
% for t=1:step:length(tlist)
% disp(sprintf('time = %d',(t-1)*dt));
% plot(xlist,sol(t,:,1),xlist,sol(t,:,2),xlist,sol(t,:,3));
% legend('c','f','w');
% ylim([-0.1,lambda/mu]);
% % ylim([-0.1,1]);
% drawnow;
% %     pause(0.0025);
% end
% % -------------------------- %

%%% Heatmaps of solutions
c = sol(1:10:end,:,1);
f = sol(1:10:end,:,2);
w = sol(1:10:end,:,3);

[X,Y] = meshgrid(tlist(1:10:end),xlist);

% subplot(3,1,1)
% contourf(c')
% 
% subplot(3,1,2)
% contourf(f')
% 
% subplot(3,1,3)
% contourf(w')

%%% plot climax
c = round(c,8); % round small entries to 0

figure()
contourf(X',Y',c)
xlabel('Time steps','Fontsize',18)
ylabel('Space','Fontsize',18)
colorbar
shading interp
title('Climax','Fontsize',16)

%%% plot attack
figure()
contourf(X',Y',f)
xlabel('Time steps','Fontsize',18)
ylabel('Space','Fontsize',18)
colorbar
shading interp
title('Attack','Fontsize',16)

%%% plot oxygen
figure()
contourf(X',Y',w)
xlabel('Time steps','Fontsize',18)
ylabel('Space','Fontsize',18)
colorbar
shading interp
title('Oxygen','Fontsize',16)





%%% Functions =============================================================



function [D,s] = my_pde(x,t,U)
global ru rv rw bu bv su sw ss Dc Df Dw

%%% Diffusions coefficients
D = [Dc; Df; Dw];

u = U(1,:);
v = U(2,:);
w = U(3,:);

% dcdt = (r1*w.^n./(b^n + w.^n)).*c.*(1 - (c + f)) - d1*c;
% dfdt = (r2*(1 - w.^n./(b^n + w.^n))).*f.*(1 - (f + c)) - d2.*f - q.*f.*w;
% dwdt = lambda - mu*w - eta*c.*w;

dudt = -ru*u + (ss*(u.^2 + bu))./(v.*(1+su*u.^2).*(1+sw*w));
dvdt = -rv*v + ss*u.^2 + bv;
dwdt = -rw*w + rw*u;

s = [dudt; dvdt; dwdt];

end