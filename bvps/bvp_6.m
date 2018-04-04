%
% bvp4.m 
% second order finite difference method for the bvp
%   u''(x) = f(x),   u'(ax)=sigma,   u(bx)=beta
% fourth order finite difference method for the bvp
%   u'' = f,   u'(ax)=sigma,   u(bx)=beta
% Using 5-pt differences on an arbitrary grid.
% Should be 4th order accurate if grid points vary smoothly.
%
% Different BCs can be specified by changing the first and/or last rows of 
% A and F.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter2  (2007)

clc;
clear all;
close all;

ax = 0;
bx = 3;
sigma = -5;   % Dirichlet boundary condition at ax
beta = 3;     % Neumann boundary condtion at bx

f = @(x) exp(x);  % right hand side function
utrue = @(x) exp(x) + (beta-exp(bx))*(x - ax) + sigma - exp(ax);  % true soln

% true solution on fine grid for plotting:
xfine = linspace(ax,bx,101);
ufine = utrue(xfine);

% Solve the problem for ntest different grid sizes to test convergence:
%m1vals = [10 20 40 80];
%m1vals = [40 80 100 120];
m1vals = [80 100 120 140 160 180 200];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);      % to hold errors

fprintf('\nChoose the finite difference method approx. right boundary condition :\n');
FDMethod = input('1) One-sided 4th order \n2) One-sided 6th order \n\n');

for jtest=1:ntest
  m1 = m1vals(jtest);
  m2 = m1 + 1;
  m = m1 - 1;                 % number of interior grid points
  hvals(jtest) = (bx-ax)/m1;  % average grid spacing, for convergence tests

  % set grid points:  
  gridchoice = 'uniform';
  x = xgrid(ax,bx,m,gridchoice);   

  % set up matrix A (using sparse matrix storage):
  A = spalloc(m2,m2,7*m2);   % initialize to zero matrix

  % first row for Neumann BC on u'(x(1))
  A(1,1:2) = fdcoeffF(0, x(1), x(1:2));
  % second row for u''(x(2))
  A(2,1:8) = fdcoeffF(2, x(2), x(1:8));
  % third row for u''(x(3))
  A(3,1:8) = fdcoeffF(2, x(3), x(1:8));

  % interior rows:
  for i=4:m-1
     A(i,i-3:i+3) = fdcoeffF(2, x(i), x((i-3):(i+3)));
  end

  % row m for u''(m)
  A(m,m-5:m2) = fdcoeffF(2,x(m),x(m-5:m2));
  % next to last row for u''(x(m+1))
  A(m1,m-5:m2) = fdcoeffF(2,x(m1),x(m-5:m2));
  % last row for Neumann BC on u(x(m+2))

  if FDMethod == 1
      A(m2,m-2:m2) = fdcoeffF(1,x(m2),x(m-2:m2));
  elseif FDMethod == 2
      A(m2,m-4:m2) = fdcoeffF(1,x(m2),x(m-4:m2));
  end

  % Right hand side:
  F = f(x); 
  F(1) = sigma;  
  F(m2) = beta;
  
  % solve linear system:
  U = A\F;


  % compute error at grid points:
  uhat = utrue(x);
  err = U - uhat;
  E(jtest) = max(abs(err));  
  disp(' ')
  disp(sprintf('Error with %i points is %9.5e',m2,E(jtest)))

  figure
  clf
  plot(x,U,'o')  % plot computed solution
  title(sprintf('Computed solution with %i grid points',m2));
  hold on
  plot(xfine,ufine)  % plot true solution
  hold off
  
  % pause to see this plot:  
  drawnow
  % input('Hit <return> for next plot ');
  
  end

error_table(hvals, E);   % print tables of errors and ratios
figure
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit
