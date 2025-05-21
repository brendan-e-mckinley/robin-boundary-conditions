set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clear;clc

size = 1;
n = 4; 
N = 10*2.^n;

uApproxArray = zeros(N(size)^2, size);
qArray = zeros(N(size)-1, size);
uExactArray = zeros(N(size)^2, size);

for j = 1:size
    Nsquared = N(j)^2;
    x = linspace(0,2*pi,N(j)+2)';
    x_int = x(2:end-1);
    dx = x_int(2)-x_int(1);
    dib = dx;
    [X,Y] = meshgrid(x_int,x_int);
    f = zeros(N(j));
    
    % Gaussian approximation of delta function
    delta_a = @(r,a) (1/(2*pi*a^2))*exp(-0.5*(r/a).^2);
    delta = @(r) delta_a(r,1.2*dx);
    
    % Immersed boundary
    theta = linspace(0,2*pi,N(j)+1)';
    theta_int = theta(1:end-1);
    xib = pi + cos(theta_int);
    yib = pi + sin(theta_int);
    rhs2 = 2*dib*cos(2.*xib).*cos(2.*yib);

    % construct D2 matrix (if robin conditions at the boundaries)
    %main = -2.*ones(N(j), 1);
    %main(1) = 1/(2*dx) - 2;
    %main(N(j)) = -1/(2*dx) - 2;
    %upper = ones(N(j), 1);
    %lower = ones(N(j), 1);
    %D2 = (1/dx^2)*spdiags([lower, main, upper], -1:1, N(j), N(j));
    %D2(1, N(j)) = 1/(2*dx^3);
    %D2(N(j), 1) = -1/(2*dx^3);
    %In = speye(N(j));
    %Lap = kron(In,D2) + kron(D2,In);

    RHS = [f(:); rhs2(:)];
    
    tol = .01*dx*dx;
    D = @(u) applyLap(u,dx,N(j));
    S = @(q) spreadQ(X,Y,xib,yib,q,delta);
    J = @(u) assembleJ(X,Y,theta_int,u,delta,dib);
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

uib = interpPhi(X,Y,xib,yib,u_approx,delta);

figure(1);
numericalSurf = surf(X,Y,u_approx);
hold on;
ib = plot3(xib, yib, uib, 'o');
set(numericalSurf, {'DisplayName'}, {['$u(x, y)$']});
set(ib, {'DisplayName'}, {['Immersed Boundary']});
set(ib, {'Color'}, {'black'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('potential', 'Interpreter', 'latex', 'FontSize', 22);
title('$u(x, y)$, numerical solution', 'Interpreter', 'latex', 'FontSize', 22);
leg = legend('FontSize', 22, 'Location', 'sw');
set(leg, 'Interpreter','latex');

% Solving for normal derivative 
normal = interpPhi(X,Y,xib,yib,u_approx,delta) - cos(2.*xib).*cos(2.*yib);
for k = 1:length(uib)
    if (mod(k, 2) == 0)
        run = 1/sqrt(1 + normal(k)^2);
        if k == 2
            normalPlot = quiver3(xib(k), yib(k), uib(k), run*cos(theta(k)), run*sin(theta(k)), normal(k)/sqrt(1 + normal(k)^2), 'Color', 'r', 'LineWidth', 2);
            set(normalPlot, {'DisplayName'}, {'$\left.\frac{\partial u}{\partial n}\right|_{\Omega}$ '});
        else
            quiver3(xib(k), yib(k), uib(k), run*cos(theta(k)), run*sin(theta(k)), normal(k)/sqrt(1 + normal(k)^2), 'Color', 'r', 'LineWidth', 2, 'HandleVisibility','off');
        end
    end
    hold on
end

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

function Jfull = assembleJ(X,Y,theta_int,Phi,delta, dib)
    xib = pi + cos(theta_int);
    yib = pi + sin(theta_int);
    xouter = pi + (1 + dib)*cos(theta_int);
    youter = pi + (1 + dib)*sin(theta_int);
    xinner = pi + (1 - dib)*cos(theta_int);
    yinner = pi + (1 - dib)*sin(theta_int);
    Jfull = 2*dib*interpPhi(X,Y,xib,yib,Phi,delta) + interpPhi(X,Y,xinner,yinner,Phi,delta) - interpPhi(X,Y,xouter,youter,Phi,delta); 
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
