

% The following code serves to produce animated graphs for the 
% normal modes of vibration of a vibrating membrane.

%  Data
    L=30;
    Nt=600;             % Maximum # of time iterations.
    Nx=100;             % Number of points in the x-direction
    deltat=0.1;
    deltax=L/Nx;
    M=30
  
% Generating the grid in the xy-plane. Two arrays.
    for i=1:Nx+1
        x(i)= (i-1)*deltax;
    end
    
 % Loop in time where the surfaces are created
    for nt=1:Nt
        t=(nt-1)*deltat;
        for i=1:Nx+1
            u(i)=0;
            for n=1:M 
                u(i)=u(i)+9/(n^2*pi^2)*sin(n*pi/3)*sin(n*pi/30*x(i))*cos(2*n*pi/30*t);
            end
        end
        plot(x,u);
        axis([0 30 -1 1]);
        xlabel('x-axis');
%       view(-63,10);               % You can change these view by rotating your surface and identifying new locations.
%        pause
%       pause(0.1);
        pause(0.00001)  % The animation is produced by introducing this pause. It's that simple!
    end