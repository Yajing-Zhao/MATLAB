function [x] = TRIDIAG(A,b)
[m,n] = size(A);
L(1,1) = A(1,1);
U(n,n) = 1;
x = zeros(n,1);
y = zeros(n,1);

for i=1:n-1
    U(i,i) = 1;
    U(i,i+1) = A(i,i+1)/L(i,i);
    L(i+1,i) = A(i+1,i);
    L(i+1,i+1) = A(i+1,i+1)-L(i+1,i)*U(i,i+1);
end
y(1) = b(1)/L(1,1);
for i=2:n
    y(i) = (b(i)-L(i,i-1)*y(i-1))/L(i,i);
end


x(n) = y(n);
for i=n-1:-1:1
    x(i) = y(i) - U(i,i+1)*x(i+1);
end



