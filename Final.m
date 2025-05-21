set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clear;clc

size = 1;
n = 2; 
N = 10*2.^n;

function u = u(x, y, dx)
    u = cos(2.*x).*sin(2.*y)./(dx^2);
end

dt = 0.01;
t = .2;
NTimeSteps = t/dt;

uApproxArray = zeros(N(size)^2, size);
qArray = zeros(N(size)-1, size);
uExactArray = zeros(N(size)^2, size);

Nsquared = N(1)^2;
x = linspace(0,2*pi,N(1)+2)';
x_int = x(2:end-1);
dx = x_int(2)-x_int(1);
[X,Y] = meshgrid(x_int,x_int);
f = zeros(N(1));
uMatrix = reshape(uApproxArray(1:Nsquared, 1), 40, 40);
figure(1);
numericalSurf = surf(X,Y, uMatrix);
for k = 1:NTimeSteps
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
    
    % Immersed boundary
    theta = linspace(0,2*pi,N(j)+1)';
    theta_int = theta(1:end-1);
    xib = k/10*pi + cos(theta_int);
    yib = k/10*pi + sin(theta_int);
    uib = cos(2.*xib).*sin(2.*yib);

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
uMatrix = reshape(uApproxArray(1:Nsquared, 1), 40, 40);
set(numericalSurf, 'ZData', uMatrix);
drawnow;
end

xlim([0, 2*pi]);
ylim([0, 2*pi]);
u_approx = reshape(uApproxArray(1:N(size)^2, size),N(size),N(size));
x = linspace(0,2*pi,N(size)+2)';
x_int = x(2:end-1);
[X,Y] = meshgrid(x_int,x_int);
theta = linspace(0,2*pi,N(j)+1)';
theta = theta(1:end-1);
xib = pi + cos(theta);
yib = pi + sin(theta);
uib = cos(2.*xib).*sin(2.*yib);
figure(2);
numericalSurf = surf(X,Y,u_approx);
hold on;
ib = plot3(xib, yib, uib, 'o');
set(numericalSurf, {'DisplayName'}, {['$\Phi$']});
set(ib, {'DisplayName'}, {['Immersed Boundary']})
set(ib, {'Color'}, {'black'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('potential', 'Interpreter', 'latex', 'FontSize', 22);
title('$u(x, y)$, numerical solution', 'Interpreter', 'latex', 'FontSize', 22);
leg = legend('FontSize', 22, 'Location', 'sw');
set(leg, 'Interpreter','latex');
hold off

figure(2);
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
