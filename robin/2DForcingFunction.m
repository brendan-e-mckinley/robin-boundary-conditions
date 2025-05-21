set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clear;clc

size = 7;
n = 1:7; 
N = 10*2.^n;

uApproxArray = zeros(N(size)^2, size);
uExactArray = zeros(N(size)^2, size);

function u = u(x, y, dx)
    u = sin(x)*sin(y)/(dx^2);
end

for j = 1:size
    Nsquared = N(j)^2;
    x = linspace(0,2*pi,N(j)+2)';
    x_int = x(2:end-1);
    dx = x_int(2)-x_int(1);
    
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
    [X,Y] = meshgrid(x_int,x_int);
    f = -2.*sin(X).*sin(Y);
    uExact = sin(X).*sin(Y);
    uExactArray(1:Nsquared, j) = uExact(:);
    
    % BCs = zeros(N(j));
    % BCs(1,:) = BCs(1,:) - u(x_int, 0, dx)';
    % BCs(:,1) = BCs(:,1) - u(0, x_int, dx);
    % BCs(:,N(j)) = BCs(:,N(j)) - u(5, x_int, dx);
    % BCs(N(j),:) = BCs(N(j),:) - u(x_int, 5, dx)';
    
    f_full = f;
    sol = Lap\f_full(:);
    tol = 0.1*dx*dx; %try = 1e-14;
    %[u_sparse_gmres,FLAG,RELRES,ITER] = gmres(Lap,f(:),[],tol);
    %disp(['GMRES converged in ' num2str(ITER(2)) ' itterations and took '])
    
    %figure(1);
    %surf(X,Y,sol);
    
    %uApproxArray(1:Nsquared, j) = u_sparse_gmres;
    uApproxArray(1:Nsquared, j) = sol;
    %Rel_L2_err = sqrt(mean((u_exact(:) - u_approx(:)).^2))./sqrt(mean(u_exact(:).^2));
    %disp(['sparse relative error: ' num2str(Rel_L2_err)]);
end

u_approx = reshape(uApproxArray(1:N(3)^2, 3),N(3),N(3));
x = linspace(0,5,N(3)+2)';
x_int = x(2:end-1);
[X,Y] = meshgrid(x_int,x_int);
figure(1);
numericalSurf = surf(X,Y,u_approx);
set(numericalSurf, {'DisplayName'}, {'u_n(x, y)'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('$u(x,y)$', 'Interpreter', 'latex', 'FontSize', 22);
title('$u_n(x, y)$, numerical solution', 'Interpreter', 'latex', 'FontSize', 22);
legend('FontSize', 22, 'Location', 'sw');

uExact = reshape(uExactArray(1:N(3)^2, 3),N(3),N(3));
figure(2);
analyticSurf = surf(X,Y,uExact);
set(analyticSurf, {'DisplayName'}, {'u_a(x, y)'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('$u(x,y)$', 'Interpreter', 'latex', 'FontSize', 22);
title('$u_a(x, y)$, analytic solution', 'Interpreter', 'latex', 'FontSize', 22);
legend('FontSize', 22, 'Location', 'sw');

% calculate error
err2Norm = zeros(1, size);
for i = 1:size
    err2Norm(i) = norm(uExactArray(1:N(i)^2,i) - uApproxArray(1:N(i)^2,i)) / norm(uExactArray(1:N(i)^2,i));
end

x = [0.64 0.64];
y = [0.7 0.52];

slopeCheck2 = polyfit(log10(N), log10(err2Norm), 1);
slope2String = 'slope $\approx$ ' + string(round(slopeCheck2(1)));
polyfitString = '$\log_{10}{(Err)} = ' + string(slopeCheck2(1)) + '\log_{10}{(N)} + ' + string(slopeCheck2(2)) + '$';

figure(3);
hold off;
plotErrInf = loglog(N, err2Norm, '--o');
set(plotErrInf, {'DisplayName'}, {'Err'})
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Err$', 'Interpreter', 'latex', 'FontSize', 22);
title('Relative error of $u_n(x, y)$, $p = 2$ norm', 'Interpreter', 'latex', 'FontSize', 22);
legend('FontSize', 22, 'Location', 'ne');
annotation("textarrow", x, y, 'String', slope2String, 'Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', polyfitString, 'Interpreter', 'latex');
saveas(gcf, 'errorPoisson.png');
