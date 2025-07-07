% PS5- M matrix and plots

clear all; clc;

a      = 0.2;
phi    = 1;
delta  = 0.05;
rho    = 0.01;
sigma_ss = 0.2;
b      = 0.05;
nu     = 0.02;
sigmaN = 0.05;      % aggregate vol. of numéraire
xi     = 1;
theta  = 0.5;

N         = 200;
sigma_min = 0.01;
sigma_max = 0.5;
sigma_grid= linspace(sigma_min, sigma_max, N)';

mu_sigma  = b*(sigma_ss - sigma_grid);
sig_sigma = nu*sqrt(sigma_grid);

M = build_M(sigma_grid, mu_sigma, sig_sigma);

% Value–Function Iteration (Implicit)
% Initial guess
vartheta = 0.5*ones(N,1);

% Drift term u(vartheta)
u_fun = @(v,s) -v.*(1-v).*sigmaN.*xi.*(1-theta).*s ...
              + 0.5*v.*(1-v).*(1-2*v).*(xi.*(1-theta).*s).^2;

dt       = 0.01;
max_iter = 1e4;
tol      = 1e-6;

for it=1:max_iter
    u_new = u_fun(vartheta, sigma_grid);
    rhs   = dt*u_new + vartheta;
    LHS   = (eye(N)*(1+rho*dt) - dt*M);
    vartheta_next = LHS\rhs;
    if max(abs(vartheta_next - vartheta)) < tol
        fprintf('Converged in %d iterations\n', it);
        vartheta = vartheta_next;
        break;
    end
    vartheta = vartheta_next;
end

% Compute Equilibrium Objects
qK           = ones(N,1);
qB           = vartheta ./ (1 - vartheta);
rf           = rho - u_fun(vartheta,sigma_grid)./vartheta;
varsigma_agg = sigmaN * ones(N,1);
varsigma_idio= sigma_grid;

figure;

% ϑ
subplot(3,2,1);
plot(sigma_grid, vartheta,'b');
xlabel('$\tilde{\sigma}_t$','Interpreter','latex');
ylabel('$\vartheta$','Interpreter','latex');
title('Portfolio Share $\vartheta$','Interpreter','latex');

% q^B
subplot(3,2,2);
plot(sigma_grid, qB,'r');
xlabel('$\tilde{\sigma}_t$','Interpreter','latex');
ylabel('$q^B$','Interpreter','latex');
title('Bond Price $q^B$','Interpreter','latex');

% q^K
subplot(3,2,3);
plot(sigma_grid, qK,'g');
xlabel('$\tilde{\sigma}_t$','Interpreter','latex');
ylabel('$q^K$','Interpreter','latex');
title('Capital Price $q^K$','Interpreter','latex');

% r^f
subplot(3,2,4);
plot(sigma_grid, rf,'k');
xlabel('$\tilde{\sigma}_t$','Interpreter','latex');
ylabel('$r^f$','Interpreter','latex');
title('Risk-free Rate $r^f$','Interpreter','latex');

% Aggregate price of risk
subplot(3,2,5);
plot(sigma_grid, varsigma_agg,'m');
xlabel('$\tilde{\sigma}_t$','Interpreter','latex');
ylabel('$\varsigma$','Interpreter','latex');
title('Aggregate Price of Risk $\varsigma$','Interpreter','latex');

% Idiosyncratic price of risk
subplot(3,2,6);
plot(sigma_grid, varsigma_idio,'c');
xlabel('$\tilde{\sigma}_t$','Interpreter','latex');
ylabel('$\tilde{\varsigma}$','Interpreter','latex');
title('Idiosyncratic Price of Risk $\tilde{\varsigma}$','Interpreter','latex');

% Tidy up & save
set(gcf,'PaperPositionMode','auto');
print('PS5_plots','-dpdf','-r300');
