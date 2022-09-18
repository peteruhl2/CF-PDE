% Example of the periodic BC solver
% Solving the (two component) Fitzhugh-Nagumo equations
% dv/dt = v-v^3/3-w + D d^2 v / dx^2
% dw/dt = 0.08(v-0.8w+0.65)
%
% We give the initial condition a kick at
% x=1.  This leads to the formation of
% two traveling waves which eventually annihilate.

function ex2()
   xlist=0:0.025:4;
   nx=length(xlist);
   ic=[-1.2*ones(1,nx);
       -0.62*ones(1,nx)];
   
   % Give it a kick
   ic(1,round(nx/4)-1)=0;
   ic(1,round(nx/4))=1;
   ic(1,round(nx/4)+1)=0;
   
   dt=0.1;
   tlist=0:dt:130;
   
   disp('solving pde');
   sol=pbcpdeSolver(@fpde,ic,xlist,tlist);
   disp('done');
   
   % ---- Plot the solution --- %
   step=1/dt;  %step every second
   for t=1:step:length(tlist)
    disp(sprintf('time = %d',(t-1)*dt));
    plot(xlist,sol(t,:,1),xlist,sol(t,:,2));
    legend('v','w');
    ylim([-2,2]);
    drawnow;
%     pause(0.0025);
   end
   % -------------------------- %
end

% Here is where the pde is defined
function [D,s]=fpde(x,t,u)
   D=[0.001;0];  % Diffusion coefficients
   v=u(1,:);  w=u(2,:);
   dvdt=v-v.^3/3-w;
   dwdt=0.08*(v-0.8*w+0.65);
   s=[dvdt; dwdt];
end