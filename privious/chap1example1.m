
% Example 1.1 
% Table 1.1 and Figure 1.2 of  http://www.amath.washington.edu/~rjl/fdmbook

uptrue = -4*sin(2);
hvals = logspace(-1,-4,13);
Epu = [];  Emu = [];  E0u = [];  E3u = [];

% table headings:
disp(' ')
disp('       h              Dpu            Dmu             D0u             D3u')

for i=1:length(hvals)
   h = hvals(i);
   % approximations to u'(1):
   Dpu = fdstencil(2,-2:2);
   Dmu = (sin(1) - sin(1-h))/h;
   D0u = (sin(1+h) - sin(1-h))/(2*h);
   xpts = 1 + h*(-2:1)';
   D3u = fdcoeffF(1,1,xpts) * sin(xpts) ;
   % errors:
   Epu(i) = Dpu - uptrue;
   Emu(i) = Dmu - uptrue;
   E0u(i) = D0u - uptrue;
   E3u(i) = D3u - uptrue;

   % print line of table:
   disp(sprintf('%13.4e   %13.4e   %13.4e   %13.4e  %13.4e',...
                 h,Epu(i),Emu(i),E0u(i),E3u(i)))
   end

% plot errors:
clf
loglog(hvals,abs(Epu),'o-')
axis([5e-4 .2 1e-12 1])
hold on
loglog(hvals,abs(E0u),'o-')
loglog(hvals,abs(E3u),'o-')
hold off
