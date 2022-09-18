% Example of the periodic BC solver
% Here we solve the following pde
% One component pde containing a decay and diffusion term.
% du/dt = -u + d^2 u/dx^2
%
% The initial condition is the sinusoid
% u(t=0) = cos(x + pi/4)+1.
%
% We then plot the solution as a space-time plot
function sol = ex1()
   xlist=0:0.01:2*pi;
   ic=cos(xlist+pi/4)+1;
   tlist=0:0.1:4;
   
   % Here is where we call the solver.
   sol=pbcpdeSolver(@fpde,ic,xlist,tlist);

   % ---- Plot the solution --- %
   % For easier visualization, downsample the solution
   xPlot = xlist(1:20:end);
   solPlot = sol(:,1:20:end);
   surf(xPlot,tlist,solPlot);
   xlabel('x');
   ylabel('time');
   % -------------------------- %
   
end

%Here is where we define the diffusion and reaction parts.
function [D,s]=fpde(x,t,u)
   D=1;  % Diffusion coefficient
   s=-u; % constant decay
end