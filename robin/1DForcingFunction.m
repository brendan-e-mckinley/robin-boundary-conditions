set(0,'defaultLineLineWidth',2);
set(0,'defaultAxesFontSize',14);
set(0,'defaulttextInterpreter','latex');

clc;
clear;

% evaluate for different values of N 
size = 8;
n = 1:8; 
N = 10*2.^n;
uArray = zeros(N(size), size);
solArray = zeros(N(size), size);
for j = 1:size
    % solve u''(x) = f(x) on domain pi/4 <= x <= 9pi/4
    % BCs: u(pi/4) = u(9pi/4) = u'(pi/4)
    x = linspace(pi/4, 9*pi/4, N(j));
    dx = x(2)-x(1);
    % robin BCs: ignore endpoints for now
    x_interior = x(2:end-1);
    f = -sin(x_interior');

    % construct D2 matrix
    main = -2.*ones(N(j)-2, 1);
    main(1) = 1/(2*dx) - 2;
    main(N(j) - 2) = -1/(2*dx) - 2;
    upper = ones(N(j)-2, 1);
    lower = ones(N(j)-2, 1);

    D2 = (1/dx^2)*spdiags([lower, main, upper], -1:1, N(j)-2, N(j)-2);
    D2(1, N(j) - 2) = -1/(2*dx^3);
    D2(N(j) - 2, 1) = 1/(2*dx^3);

    RHS = f;
    % solve
    u = D2\RHS;
    %u = sin(pi/4) + u;
    % add endpoints of solution back in
    u_full = [cos(pi/4); u; cos(pi/4)];
    uArray(1:N(j), j) = u_full;
    % analytic solution
    solArray(1:N(j), j) = sin(x);
end
% compare
x = linspace(0, 5, N(3));
dx = x(2)-x(1);
figure(1);
analyticPlot = plot(x, solArray(1:N(3), 3), '-', 'linewidth', 4);
hold on
numericalPlot = plot(x, uArray(1:N(3), 3), '--o','linewidth',1,'markersize',10);
set(analyticPlot, {'DisplayName'}, {'Test Solution, u_t(x)'});
set(numericalPlot, {'DisplayName'}, {'Numerical Solution, u_n(x)'});
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$u(x)$', 'Interpreter', 'latex', 'FontSize', 22);
title('$u(x)$, analytic and numerical solutions', 'Interpreter', 'latex', 'FontSize', 22);
legend('FontSize', 22, 'Location', 'sw');

% calculate error
err2Norm = zeros(1, size);
errInfNorm = zeros(1, size);
for i = 1:size
    err2Norm(i) = norm(solArray(1:N(i),i) - uArray(1:N(i),i)) / norm(solArray(1:N(i),i));
end

for i = 1:size
    errInfNorm(i) = norm(solArray(1:N(i),i) - uArray(1:N(i),i), Inf) / norm(solArray(1:N(i),i), Inf);
end

slopeCheck2 = polyfit(log10(N), log10(err2Norm), 1);
slope2String = 'slope $\approx$ ' + string(round(slopeCheck2(1)));
polyfitString = '$\log_{10}{(Err)} = ' + string(slopeCheck2(1)) + '\log_{10}{(N)} + ' + string(slopeCheck2(2)) + '$';

slopeCheckInf = polyfit(log10(N), log10(errInfNorm), 1);
slopeInfString = 'slope $\approx$ ' + string(round(slopeCheckInf(1)));

x = [0.64 0.64];
y = [0.65 0.46];

figure(2);
hold off;
plotErr2 = loglog(N, err2Norm, '--o');
set(plotErr2, {'DisplayName'}, {'Err'})
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$Err$', 'Interpreter', 'latex', 'FontSize', 22);
title('Relative error of $u_n(x)$, $p = 2$ norm', 'Interpreter', 'latex', 'FontSize', 22);
legend('FontSize', 22, 'Location', 'ne');
annotation("textarrow", x, y, 'String', slope2String, 'Interpreter', 'latex', 'FontSize', 22);
annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', polyfitString, 'Interpreter', 'latex');

