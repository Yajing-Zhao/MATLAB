function of = outfunc(fh, h, j)

of = @(x) for i = 1:5
    x = fh(x) + h + j;
    
    
    % print out stencil:

disp(' ')
disp(sprintf('The derivative u^(%i) of u at x0 is approximated by',k))
disp(' ')
disp(sprintf('          1/h^%i * [',k))
for i=1:n-1
  if j(i) < 0
    disp(sprintf('                     %22.15e * u(x0%i*h) + ',c(i),j(i)))
  elseif j(i) == 0
    disp(sprintf('                     %22.15e * u(x0) + ',c(i)))
  else
    disp(sprintf('                     %22.15e * u(x0+%i*h) + ',c(i),j(i)))
  end
end

disp(sprintf('                     %22.15e * u(x0+%i*h) ]   ',c(n),j(n)))


% determine dominant terms in truncation error and print out:

err0 = c*(j(:).^n) / factorial(n);
err1 = c*(j(:).^(n+1)) / factorial(n+1);
if (abs(err0)) < 1e-14,  err0 = 0; end   % for centered approximations, expect
if (abs(err1)) < 1e-14,  err1 = 0; end   % one of these to be exactly 0.
disp(' ')
disp('For smooth u,')
disp(sprintf('       Error = %g * h^%i*u^(%i) + %g * h^%i*u^(%i) + ...',err0,n-k,n,err1,n-k+1,n+1))

disp(' ')