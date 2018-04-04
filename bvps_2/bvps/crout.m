function [L,U] = crout(A,n)  % returns two matrices

for i = 1:n
    L(i,1) = A(i,1);
end

for j = 1:n
    U(1,j) = A(1,j)/L(1,1);
end

for j = 2:n
    for i = j:n
        sum = 0.0;
        for k = 1:(j-1)
            sum = sum + L(i,k) * U(k,j);
        end
        L(i,j) = A(i,j) - sum;
    end

    U(j,j) = 1;

    for i = (j+1):n
        sum = 0.0;
        for k = 1:(j-1)
            sum = sum + L(j,k) * U(k,i);
        end
        U(j,i) = (A(j,i) - sum)/L(j,j);
    end
end