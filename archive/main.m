k = 2;
j = -2:2;
xbar = 1;

fh = fdextension(k,j, xbar, @(x)(sin(2*x)));

uptrue = -4*sin(2);
hvals = logspace(-1,-4,13);
Error = []; 

% table headings:
disp(' ')
disp('       h              D0u')

for i=1:length(hvals)
   h = hvals(i);
   % approximations to u'(1):
   D = fh(h);
   % errors:
   Error(i) = D - uptrue;
  
   % print line of table:
   disp(sprintf('%13.4e   %13.4e',...
                 h,Error(i)));
end

error_loglog(hvals,Error)

