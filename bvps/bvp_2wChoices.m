%
% bvp_2.m 
% second order finite difference method for the bvp
%   u''(x) = f(x),   u'(ax)=sigma,   u(bx)=beta
% Using 3-pt differences on an arbitrary nonuniform grid at the interior
% points. But, you can choose between first and second order one-sided 
% approximation at the left boundary.

% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)
% Modified by Vianey Villamizar (2017)

clc;
close all;
clear all;

ax = 0;
bx = 3;
sigma = -5;   % Dirichlet boundary condition at ax
beta = 3;     % Neumann boundary condtion at bx

f = @(x) exp(x);  % right hand side function
utrue = @(x) exp(x) + (beta-exp(bx))*(x - ax) + sigma - exp(ax);  % true soln

% true solution on fine grid for plotting:
xfine = linspace(ax,bx,101);
ufine = utrue(xfine);

fprintf('\nChoose the finite difference method approx. right boundary condition :\n');
FDMethod = input('1) One-sided 2nd order \n2) Centered 2nd order \n\n');

% Solve the problem for ntest different grid sizes to test convergence:
%m1vals = [10 20 40 80];
m1vals = [40 80 120 160];
%m1vals = [5 10 20 40];
ntest = length(m1vals);
hvals = zeros(ntest,1);  % to hold h values
E = zeros(ntest,1);      % to hold errors

for jtest=1:ntest
  m1 = m1vals(jtest);
  m2 = m1 + 1;
  m = m1 - 1;                 % number of interior grid points
  hvals(jtest) = (bx-ax)/m1;  % average grid spacing, for convergence tests

  % set grid points:  
  gridchoice = 'uniform';          % see xgrid.m for other choices
  x = xgrid(ax,bx,m,gridchoice);   

   % Forcing Term
  F = f(x); 
  F(1) = sigma;
  
  % set up matrix A (using sparse matrix storage):
  A = spalloc(m2,m2,3*m2);   % initialize to zero matrix

  % Construction of first row of matrix A
  if FDMethod == 1
        % Last row for One-Sided second Order approx. for Neumann BC
        % at bx (use 3 points):
        A(m2,m:m2) = fdcoeffF(1,x(m2),x(m:m2));
        F(m2) = beta;
        
  elseif FDMethod == 2
        % Last row for Centered Second Order approx. for Neumann BC 
        % at bx (use 3 points):
        A(m2,m1:m2) = fdcoeffF(1,x(m2),x(m1:m2));
        F(m2) = beta - f(bx)*hvals(jtest)/2;
  end
  
 % interior rows:
  for i=2:m1
      A(i,i-1:i+1) = fdcoeffF(2, x(i), x((i-1):(i+1)));
  end

 % last row for Dirichlet BC at bx:
  % A(m2,m:m2) = fdcoeffF(0,x(m2),x(m:m2)); 
  A(1,1:3) = fdcoeffF(0,x(1),x(1:3));
        

  
   
  FA = full(A);
  % solve linear system:
  U = A\F;

  % compute error at grid points:
  uhat = utrue(x);
  err = U - uhat;
  E(jtest) = max(abs(err));  
  disp(' ')
  disp(sprintf('Error with %i points is %9.5e',m2,E(jtest)));

  figure;
  plot(x,U,'o')  % plot computed solution
  title(sprintf('Computed solution with %i grid points',m2));
  hold on
  plot(xfine,ufine)  % plot true solution
  %hold off
  
  % pause to see this plot:  
  drawnow
  %input('Hit <return> to continue ');
end

figure;
error_table(hvals, E);   % print tables of errors and ratios
error_loglog(hvals, E);  % produce log-log plot of errors and least squares fit

