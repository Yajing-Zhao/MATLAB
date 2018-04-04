function error_loglog(h,E)
%
% Produce log-log plot of E vs. h.
% Estimate order of accuracy by doing a linear least squares fit.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/  (2007)

h = h(:);            % make sure it's a column vector
E = E(:);            % make sure it's a column vector
ntest = length(h);
clf
loglog(h,E,'o-')
axis([.5*min(h) 1.5*max(h)  .5*min(E) 1.5*max(E)])
title('log-log plot of errors vs. h')
hold on