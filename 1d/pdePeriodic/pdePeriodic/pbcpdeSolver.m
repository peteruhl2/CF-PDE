function sol = pbcpdeSolver(fpde,ic,xlist,tlist)
% pbcpdeSolver - Periodic boundary condition reaction diffusion solver
% Usage:
%    sol = pbcpdeSolver(fpde,ic,xlist,tlist);
% Returns:
%    sol is a three dimensional array (nt X nx+1 X npde),
%    where nt, nx+1, and npde are respectively the number of time points,
%    position points, and number of reaction diffusion equations.
% Arguments:
%    fpde is a function of the form
%       [D,s] = fpde(x,t,u);
%          x can be a row vector (1 X nx) of grid points
%          t is a scalar time variable
%          u is a matrix (npde X nx) containing the values of the RD eqn.
%       and
%          D is a column vector (npde X 1) for the diffusion coefficients.
%          s is a matrix (npde X nx) representing the source term
%    ic are the initial conditions (npde X nx)
%    xlist is a vector (1 x nx+1) of monotonicaly increasing values where
%       the concentration at xlist(nx+1) is same as at xlist(1).
%    tlist is the list of times for which to store the solution.
%
% Written by Albert J. Bae
% Last modified: March 19, 2014

sol=[];
if nargin ~= 4
  error('pdeSolver requires 4 inputs')
end

t = tlist(:);
nt = length(t);
if nt < 3
  error('Not enough time points')
end
if any(diff(t) <= 0)
  error('tlist must be increasing')
end

x = xlist(:);
nx = length(x);  % Note that nx in the code here is nx+1 in
                 % the description above.
if nx < 4
  error('not enough points in the grid')
end
h = diff(x);
if any(h <= 0)
  error('xlist must be increasing')
end

sz=size(ic);
if(sz(2)~=nx)
    error('initial condition has incorrect dimensions');
end
npde=sz(1);

% Check the output of pde function
[D, S] = fpde(x,t(1),ic);
sz2=size(D); sz3=size(S);
if(sz2(1)~=npde||sz2(2)~=1)
    error('RD diffusion term should be a column vector');
end
if(sz3(1)~=sz(1)||sz3(2)~=sz(2))
    error('RD source term has incorrect dimensions.');
end


% Remove the last (periodic) point and
% collapse the matrix down to a row vector
n=nx-1;
ic2 = ic(:,1:n);
ic2 = reshape(ic2',1,[]);

x2  = repmat(x(1:n)',[1,npde]);
h  = repmat(h',[1,npde]);
ind = repmat(1:n,[npde,1]);
for i=1:npde-1
    ind(i+1,:)=ind(i+1,:)+i*n;
end
iL=circshift(ind,[0,1]);
iR=circshift(ind,[0,-1]);
iL=reshape(iL',1,[]);
iR=reshape(iR',1,[]);

s  = ode15s(@pdeodes,t',ic2);
sz = size(s.y);
sol=zeros([length(t),n+1,npde]);

for i=1:npde
    for j=1:n
        k=(i-1)*n+j;
        sol(:,j,i)=interp1(s.x,s.y(k,:),t);
    end
    sol(:,n+1,i)=sol(:,1,i);
end

% nested function
function dydt = pdeodes(tnow,y)
    % Get the diffusion coefficients and
    % source terms
    y=y(:)';
    y2=reshape(y,n,npde)'; % Matrixize it
    [D, S] = fpde(x(1:n),tnow,y2);
    
    % source term part
    dydt = reshape(S',1,[]);
    
    % diffusion part
    D = repmat(D,[1,n]);
    D = reshape(D',1,[]);
    Delta2 = y(iR)-y;
    Delta1 = y - y(iL);
    h1 = h(iL);
    h2 = h;
    dydt = dydt + ...
        (2*D).*(h1.*Delta2-h2.*Delta1)./(h1.*h2.*(h1+h2));
    dydt = dydt(:);
end

end
