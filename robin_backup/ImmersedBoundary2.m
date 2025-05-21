set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clear;clc

size = 1;
n = 4; 
N = 10*2.^n;

function u = u(x, y, dx)
    u = sin(2.*x).*sin(2.*y)./(dx^2);
end

uApproxArray = zeros(N(size)^2, size);
qArray = zeros(N(size)-1, size);
uExactArray = zeros(N(size)^2, size);

for j = 1:size
    Nsquared = N(j)^2;
    x = linspace(0,2*pi,N(j)+2)';
    x_int = x(2:end-1);
    dx = x_int(2)-x_int(1);
    [X,Y] = meshgrid(x_int,x_int);
    f = zeros(N(j));
    
    % Gaussian approximation of delta function
    delta_a = @(r,a) (1/(2*pi*a^2))*exp(-0.5*(r/a).^2);
    delta = @(r) delta_a(r,1.2*dx);
    
    % Normal derivative
    theta = linspace(0,2*pi,N(j)+1)';
    theta_int = theta(1:end-1);
    xib2 = pi + (1+dx)*cos(theta_int);
    yib2 = pi + (1+dx)*sin(theta_int);
    uib2 = sin(2.*xib2).*sin(2.*yib2);
    xib0 = pi + (1-dx)*cos(theta_int);
    yib0 = pi + (1-dx)*sin(theta_int);
    uib0 = sin(2.*xib0).*sin(2.*yib0);

    % Immersed boundary
    xib = pi + cos(theta_int);
    yib = pi + sin(theta_int);
    uib = (uib2 - uib0)/(2*dx);

    % construct D2 matrix
    main = -2.*ones(N(j), 1);
    main(1) = 1/(2*dx) - 2;
    main(N(j)) = -1/(2*dx) - 2;
    upper = ones(N(j), 1);
    lower = ones(N(j), 1);
    D2 = (1/dx^2)*spdiags([lower, main, upper], -1:1, N(j), N(j));
    D2(1, N(j)) = 1/(2*dx^3);
    D2(N(j), 1) = -1/(2*dx^3);
    In = speye(N(j));
    Lap = kron(In,D2) + kron(D2,In);

    RHS = [f(:); uib(:)];
    
    tol = .01*dx*dx;
    D = @(u) applyLap(u,dx,N(j));
    S = @(q) spreadQ(X,Y,xib,yib,q,delta);
    J = @(u) interpPhi(X,Y,xib,yib,u,delta);
    Vec = @(v) applyFull(v,Nsquared,N(j),D,S,J);
    [full_gmres,FLAG,RELRES,ITER] = gmres(Vec,RHS,1000,tol,1);

    uApproxArray(1:Nsquared, j) = full_gmres(1:Nsquared);
    qArray(1:N(j), j) = full_gmres(Nsquared+1:end);
end

% Only plotting one of these, so we need to recompute some parameters
xlim([0, 2*pi]);
ylim([0, 2*pi]);
u_approx = reshape(uApproxArray(1:N(size)^2, size),N(size),N(size));
x = linspace(0,2*pi,N(size)+2)';
x_int = x(2:end-1);
[X,Y] = meshgrid(x_int,x_int);

% Normal derivative
theta = linspace(0,2*pi,N(j)+1)';
theta_int = theta(1:end-1);
xib2 = pi + (1+dx)*cos(theta_int);
yib2 = pi + (1+dx)*sin(theta_int);
uib2 = sin(2.*xib2).*sin(2.*yib2);
xib0 = pi + (1-dx)*cos(theta_int);
yib0 = pi + (1-dx)*sin(theta_int);
uib0 = sin(2.*xib0).*sin(2.*yib0);

% Immersed boundary
xib = pi + cos(theta_int);
yib = pi + sin(theta_int);
uib = (uib2 - uib0)/(2*dx);

figure(1);
numericalSurf = surf(X,Y,u_approx);
hold on;
ib = plot3(xib, yib, uib, 'o');
set(numericalSurf, {'DisplayName'}, {['$\Phi$']});
set(ib, {'DisplayName'}, {['Immersed Boundary']});
set(ib, {'Color'}, {'black'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('potential', 'Interpreter', 'latex', 'FontSize', 22);
title('$u(x, y)$, numerical solution', 'Interpreter', 'latex', 'FontSize', 22);
leg = legend('FontSize', 22, 'Location', 'sw');
set(leg, 'Interpreter','latex');
%hold off

%figure(2);
xibunit = pi + 2*cos(theta_int);
yibunit = pi + 2*sin(theta_int);
for k = 1:length(uib)
    ibmag = 1/sqrt(1 + uib(k)^2);
    normal = quiver3(xib(k), yib(k), uib(k), ibmag*cos(theta(k)), ibmag*sin(theta(k)), uib(k) + uib(k)/sqrt(1 + uib(k)^2), 'Color', 'r');
    hold on
end
%ibmag = 1/sqrt(1 + uib.^2);
%U = ibmag'.*cos(theta_int);
%V = ibmag'.*sin(theta_int);
%W = uib + uib./sqrt(1 + uib.^2);

%normal = quiver3(xib, yib, uib, U, V, W, 'Color', 'r');
hold on
set(normal, {'DisplayName'}, {'Normal Derivatives'});

figure(3);
charges = plot3(xib,yib,qArray(1:N(size),size),'-o');
set(charges, {'DisplayName'}, {'$q$'});
set(charges, {'Color'}, {'black'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('$q$', 'Interpreter', 'latex', 'FontSize', 22);
title('Charge strength, $q$', 'Interpreter', 'latex', 'FontSize', 22);
leg = legend('FontSize', 22, 'Location', 'ne');
set(leg, 'Interpreter','latex');

%%
function Sq = spreadQ(X,Y,xq,yq,q,delta)
    Sq = 0*X;
    Nq = length(q);
    dx = X(1,2)-X(1,1);
    dy = Y(2,1)-Y(1,1);
    for k = 1:Nq
        Rk = sqrt((X-xq(k)).^2 + (Y-yq(k)).^2);
        Sq = Sq + q(k)*delta(Rk);
    end
    Sq = Sq(:);
end

function Jphi = interpPhi(X,Y,xq,yq,Phi,delta)
    Jphi = 0*xq;
    Phi = reshape(Phi,length(X(1,:)),length(Y(1,:))); % reshape from vector to scalar
    Nq = length(xq);
    dx = X(1,2)-X(1,1);
    dy = Y(2,1)-Y(1,1);
    for k = 1:Nq
        Rk = sqrt((X-xq(k)).^2 + (Y-yq(k)).^2);
        Jphi(k) = dx*dy*sum(sum(Phi.*delta(Rk)));
    end
end

function y = applyLap(u,dx,N)
    u = reshape(u,N,N);
    y = 0*u;
    for i = 1:N
        for j = 1:N
            u_E = 0;
            u_W = 0;
            u_N = 0;
            u_S = 0;
            u_C = u(i,j);
            if i > 1; u_W = u(i-1,j); end
            if i < N; u_E = u(i+1,j); end
            if j > 1; u_S = u(i,j-1); end
            if j < N; u_N = u(i,j+1); end
            y(i,j) = (1/(dx^2))*(u_W + u_E + u_S + u_N - 4.0*u_C);
        end
    end
    y = y(:);
end

function y = applyFull(v, N, Nib, applyLap, spreadQ, interpPhi)
    y = zeros(N + Nib, 1);
    SQ = spreadQ(v(N+1:end));
    y(1:N) = applyLap(v(1:N)) - SQ;
    y(N+1:end) = interpPhi(v(1:N));
end
