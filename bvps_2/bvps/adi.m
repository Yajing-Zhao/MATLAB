clear all;
close all;

Nx=5; %interior points in x-axis
Ny=3; %interior points in y-axis


detat=0.01; %time step
Nt = 0.05/detat;
detax=1/Nx; 
detay=1/Ny;
rx = detat/(detax)^2; %rx
ry = detat/(detay)^2; %ry
Fp = zeros((Nx-1),(Ny-1));
Fc = zeros((Nx-1),(Ny-1));
U = zeros((Nx-1),(Ny-1));

Ap =  full(gallery('tridiag',(Nx-1)*(Ny-1),-rx/2,1+rx,-rx/2)); %Ap
Ac =  full(gallery('tridiag',(Nx-1)*(Ny-1),-ry/2,1+ry,-ry/2)); %Ac

for i=1:Nx-1
        for j=1:Ny-1
            U(i,j)=0;
        end
end

% solve Ap*U(n+1/2)=Fc to get U(n+1/2)    
for n=1:Nt
    for i=1:Nx-1
        Fp(i,1)=U(i,1)+ry*(-2*U(i,1)+U(i,2))/2+detat/2;
        Fp(i,Ny-1)=U(i,Ny-1)+ry*(U(i,Ny-2)-2*U(i,Ny-1))/2+detat/2;
        for j=2:Ny-2
            Fp(i,j)=U(i,j)+ry*(U(i,j-1)-2*U(i,j)+U(i,j+1))/2+detat/2;
        end
    end
    Fp = Fp(:);
    U = TRIDIAG(Ap,Fp(:)); 
    U = reshape(U,(Nx-1),(Ny-1));
% solve Ac*U(n+1)=Fc to get U(n+1)   

    for j=1:Ny-1
        Fc(1,j)= U(1,j)+rx*(-2*U(1,j)+U(2,j))/2+detat/2;
        Fc(Nx-1,j)=U(Nx-1,j)+rx*(U(Nx-2,j)-2*U(i,j))/2+detat/2;
        for i=2:Nx-2
            Fc(i,j)=U(i,j)+rx*(U(i-1,j)-2*U(i,j)+U(i+1,j))/2+detat/2;
        end
    end
    %Fc = reshape(Fc,(Nx-1)*(Ny-1),1);
    U = TRIDIAG(Ac,Fc(:)); 
end
x = linspace(0,1,Nx-1);
y = linspace(0,1,Ny-1);

surf(x,y,U)


    
