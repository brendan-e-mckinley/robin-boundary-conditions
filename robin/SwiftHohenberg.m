%% Swift-Hohenberg

set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clear;clc

colormap('turbo');

N = 64;
Nsquared = N^2;

r = 0.7;
x = linspace(0,20*pi,N+2)';
x_int = x(2:end-1);
dx = x_int(2)-x_int(1);
epsilon = dx;
[X,Y] = meshgrid(x_int,x_int);

rng(0); % set a seed to make randomness reproducible
U = epsilon*randn(N*N,1); % random initial data with variance epsilon^2
U = U - mean(U); % make the initial data have mean = 0;

% Time grid sizes
dt = 0.01;
tFinal = 10;
Nt = tFinal / dt;
t = 0:dt:tFinal;

% laplacian matrix
main = -2.*ones(N, 1);
main(1) = 1/(2*dx) - 2;
main(N) = -1/(2*dx) - 2;
upper = ones(N, 1);
lower = ones(N, 1);

D2 = (1/dx^2)*spdiags([lower, main, upper], -1:1, N, N);
D2(1, N) = 1/(2*dx^3);
D2(N, 1) = -1/(2*dx^3);
In = speye(N);
Lap = kron(In,D2) + kron(D2,In);
InBig = speye(Nsquared);

% initial phi values
phi_current = U;
% final phi array
phi_array = zeros(Nsquared*Nt, 1);
% integral array
intArray = zeros(1,Nt);
% compute LHS 
LHS = InBig - r*dt*InBig + dt*(InBig + Lap)^2;

figure(1);
phi_current_plot = surf(X, Y, reshape(phi_current, N, N)); 
axis([0 20*pi 0 20*pi]);
set(phi_current_plot, {'DisplayName'}, {'$u(x, y, t = 0)$'})
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
title('$u$ at time $t = 0$', 'Interpreter', 'latex', 'FontSize', 22);
leg = legend('FontSize', 22, 'Location', 'ne');
set(leg, 'Interpreter','latex');
view(0,90);

for j = 1:Nt
    % numerical solution
    % compute RHS
    RHS = phi_current;
    phi_star = LHS\RHS;

    phi_current = (phi_star.*exp(dt))./sqrt(1 + phi_star.^2*(exp(2*dt) - 1));
    intArray(j) = sum(phi_current, "all")*dx^2;
    phi_array((j - 1)*Nsquared + 1:j*Nsquared) = phi_current;

    %set(phi_current_plot, 'ZData', reshape(phi_current, N, N));
    %drawnow;
end

figure(2);
tEqual50 = surf(X, Y, reshape(phi_current, N, N)); 
axis([0 20*pi 0 20*pi]);
set(tEqual50, {'DisplayName'}, {'$u(x, y, t = 10)$'})
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
title('$u$ at time $t = 10$', 'Interpreter', 'latex', 'FontSize', 22);
leg = legend('FontSize', 22, 'Location', 'ne');
set(leg, 'Interpreter','latex');
view(0,90);

% figure(3);
% mass = plot(t, intArray);
% set(mass, {'DisplayName'}, {'mass'})
% xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('$\int\phi(x, y, t) dx dy$', 'Interpreter', 'latex', 'FontSize', 22);
% title('mass over time', 'Interpreter', 'latex', 'FontSize', 22);
% leg = legend('FontSize', 22, 'Location', 'ne');
% set(leg, 'Interpreter','latex');