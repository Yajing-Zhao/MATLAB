PPW = 30;
r0 = 1;
R0 = 5;
k = 2*pi;
m = (R0-r0)*PPW;
n = int16(2*pi*r0*PPW);
I = eye(m);
A = zeros(m*n,m*n);
A(1:m,1:m) = I;
r = zeros(1,n);
theta = zeros(1,m+1);
F = zeros(n*m,1);
Uinc = exp(i*k*(r(i).*cos(theta(j))));
for i =1:m
    F(i) = -Unic % need to define F
delta_r = 0.0214;
delta_theta = 2*pi/m;

beta=zeros(n,1);
alpha=zeros(n,1);
for i=1:m+1
    theta(j) = (j-1)*delta_theta;
end
for j=1:n
    r(i) = (i-1)*delta_r;
end
for i=2:n-1
    D_minus = (2*r(i)-delta_r)/(2*r(i)*(delta_r)^2)*I;
    D_plus =(2*r(i)-delta_r)/(2*r(i)*(delta_r)^2)*I; 
    T = full(gallery('tridiag',m,beta(i),alpha(i),beta(i)));
    T(1,m) = beta(i);
    T(m,1) = beta(i);
    alpha(i) = k^2 - 2/(delta_r)^2-2/((delta_theta)^2*(r(i))^2);
    beta(i) = 1/((delta_theta)^2*(r(i))^2);
    A((i-1)*m+1:i*m,(i-2)*m+1:(i-1)*m) = (2*r(i)-delta_r)/(2*r(i)*(delta_r)^2)*I;
   
    A((i-1)*m+1:i*m,(i-1)*m+1:i*m) = full(gallery('tridiag',m,beta(i),alpha(i),beta(i)));
    A((i-1)*m+1:i*m,i*m+1:(i+1)*m)=(2*r(i)-delta_r)/(2*r(i)*(delta_r)^2)*I;     
end

A((n-1)*m+1:n*m,(n-2)*m+1:(n-1)*m)=2/((delta_r)^2)*I;

A((n-1)*m+1:n*m,(n-1)*m+1:n*m)=2/((delta_r)^2)*I;





