set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clc;
clear;

% evaluate for N = 64
N = 256;

% Time grid sizes
tSize = 5;
n = 4:8; 
Nt = 10*2.^n;
Tfinal = 1;

% solve u_t(x) = u_xx(x) on domain pi/4 <= x <= 9pi/4
% BCs: u(pi/4) = u(9pi/4) = u'(pi/4)
x = linspace(0, 2*pi, N);
dx = x(2)-x(1);
% robin BCs: ignore endpoints for now
x_interior = x(2:end-1);

% Initial conditions
u_current = sin(x_interior) + cos(x_interior);
u_current = u_current';
    
figure(1);
plotU = plot(x_interior, u_current);
axis([0 2*pi -2 2]);
for j = 1:tSize
    % time step size
    dt = Tfinal/Nt(j);
    
    % construct D2 matrix
    main = -2.*ones(N-2, 1);
    main(1) = 1/(2*dx) - 2;
    main(N - 2) = -1/(2*dx) - 2;
    upper = ones(N-2, 1);
    lower = ones(N-2, 1);

    D2 = (1/dx^2)*spdiags([lower, main, upper], -1:1, N-2, N-2);
    D2(1, N - 2) = -1/(2*dx^3);
    D2(N - 2, 1) = 1/(2*dx^3);

    for k = 1:Nt(j)
        RHS = u_current(1:N-2);
        % compute LHS 
        LHS = speye(N-2) - dt*D2;

        tol = .01*dx*dx;
        %[u_next,FLAG,RELRES,ITER] = gmres(LHS,RHS(:));
        u_next = LHS\RHS;
        u_current = u_next;
        
        set(plotU, 'YData', u_current);
        drawnow;
    end
end

figure(2);
plotU = plot(x(2:N-1), u_current);
