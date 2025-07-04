clear all;
clc;
% Parameters and grid
a_e      = 0.11;   a_h    = 0.03;   % production rates
rho_0    = 0.04;                  % time preference
rho_e_d  = 0.01; rho_h_d = 0.01;  % death rates
rho_e    = rho_0 + rho_e_d;       % expert’s discount rate
rho_h    = rho_0 + rho_h_d;       % household’s discount rate
zeta     = 0.05;                  % probability of becoming an expert
delta    = 0.05; sigma   = 0.1;   % decay rate / volatility
phi      = 10;   % adjustment cost / equity constraint

N        = 10000;                   % grid size
eta      = linspace(0.0001,0.9999,N)';  % grid for \eta

% Solution
% Solve for q (0)
q0 = (1 + a_h * phi) / (1 + rho_h * phi);

% Inner loop
[Q, SSQ, Kappa, Iota] = inner_loop_log_without_alpha(eta, q0, a_e, a_h, rho_e, rho_h, sigma, phi);

S      = (Kappa - eta) .* SSQ;        % \sigma_{\eta^e} -- arithmetic volatility of \eta^e
Sg_e   = S ./ eta;                  % \sigma^{\eta^e} -- geometric volatility of \eta^e
Sg_h   = -S ./ (1 - eta);           % \sigma^{\eta^h} -- geometric volatility of \eta^h

VarS_e = Kappa ./ eta .* SSQ;         % \varsigma^e -- experts’ price of risk
VarS_h = (1 - Kappa) ./ (1 - eta) .* SSQ;  % \varsigma^h -- households’ price of risk

CN_e   = rho_e;                     % experts’ consumption-to-networth ratio
CN_h   = rho_h;                     % households’ consumption-to-networth ratio

MU = eta .* (1 - eta) .* ((VarS_e - SSQ) .* (Sg_e + SSQ) ...
     - (VarS_h - SSQ) .* (Sg_h + SSQ) ...
     - (CN_e - CN_h) ...
     + (rho_h_d .* zeta .* (1 - eta) - rho_e_d .* (1 - zeta) .* eta) ...
       ./ (eta .* (1 - eta)));  % \mu_{\eta^e} -- arithmetic drift of \eta^e


% Create figure with specific size FOR Price and Amplification
figure('Position', [100, 100, 900, 400]);

% Find transition point where kappa becomes >= 1
transition_idx = find(Kappa >= 1, 1, 'first');
if ~isempty(transition_idx)
    transition_eta = eta(transition_idx);
else
    transition_eta = 0.4; % fallback value
end

% Plot data
subplot(1,2,1)
plot(eta, Q, 'b-', 'LineWidth', 2);
hold on
% Add vertical dashed line at transition point
plot([transition_eta transition_eta], ylim, 'k:', 'LineWidth', 1)
hold off

xlabel('\eta', 'FontSize', 14)
ylabel('q', 'FontSize', 14)
title('Price of Capital', 'FontSize', 14, 'FontWeight', 'bold')
grid off
xlim([0 1])
ylim([0.9 1.4])

% Add region labels
text(0.15, 1.22, 'Crisis Region', 'FontSize', 8, 'Color', 'red', 'Rotation', 0)
text(0.7, 1.22, 'Normal Region', 'FontSize', 8, 'Color', 'red', 'Rotation', 0)

% Add kappa labels at top
text(transition_eta/2, 1.38, '\kappa_t^e < 1', 'FontSize', 10, 'HorizontalAlignment', 'center')
text((transition_eta+1)/2, 1.38, '\kappa_t^e = 1', 'FontSize', 10, 'HorizontalAlignment', 'center')

subplot(1,2,2)
plot(eta, SSQ-sigma, 'b-', 'LineWidth', 2);
hold on
% Add horizontal dashed line at y=0
plot([0 1], [0 0], 'k--', 'LineWidth', 1)
% Add vertical dashed line at transition point
plot([transition_eta transition_eta], ylim, 'k:', 'LineWidth', 1)
hold off

xlabel('\eta', 'FontSize', 14)
ylabel('\sigma^q', 'FontSize', 14)
title('Amplification', 'FontSize', 14, 'FontWeight', 'bold')
grid off
xlim([0 1])
ylim([-0.01 0.07])

% Add region labels
text(0.15, 0.04, 'Crisis Region', 'FontSize', 8, 'Color', 'red', 'Rotation', 0)
text(0.7, 0.04, 'Normal Region', 'FontSize', 8, 'Color', 'red', 'Rotation', 0)

% Add kappa labels at top
text(transition_eta/2, 0.065, '\kappa_t^e < 1', 'FontSize', 10, 'HorizontalAlignment', 'center')
text((transition_eta+1)/2, 0.065, '\kappa_t^e = 1', 'FontSize', 10, 'HorizontalAlignment', 'center')

% Adjust subplot spacing
set(gcf, 'PaperPositionMode', 'auto')
print('price_amplification_plot', '-dpdf', '-r300')


% Create figure with specific size FOR Iota and Kappa
figure('Position', [100, 100, 800, 400]);

% Plot data
subplot(1,2,1)
plot(eta, Iota, 'b-', 'LineWidth', 2);
hold on
% Add vertical dashed line at transition point
plot([transition_eta transition_eta], ylim, 'k:', 'LineWidth', 1)
hold off

xlabel('\eta', 'FontSize', 14)
ylabel('\iota', 'FontSize', 14)
title('Investment', 'FontSize', 14, 'FontWeight', 'bold')
grid off

subplot(1,2,2)
plot(eta, Kappa, 'b-', 'LineWidth', 2);
hold on

xlabel('\eta', 'FontSize', 14)
ylabel('\kappa', 'FontSize', 14)
title('Capital share', 'FontSize', 14, 'FontWeight', 'bold')
grid off


% Adjust subplot spacing
set(gcf, 'PaperPositionMode', 'auto')
print('Investment_CapitalShare_plot', '-dpdf', '-r300')


% Create figure with specific size FOR Debt issued
figure('Position', [100, 100, 600, 400]);

% Plot data
subplot(1,1,1)
plot(eta, Kappa-eta, 'b-', 'LineWidth', 2);

xlabel('\eta', 'FontSize', 14)
ylabel('$\frac{D^{e}_{t}}{q_{t}K_{t}}$','Interpreter', 'latex','FontSize', 14)
title('Debt Issued', 'FontSize', 14, 'FontWeight', 'bold')
grid off


% Adjust subplot spacing
set(gcf, 'PaperPositionMode', 'auto')
print('Debt_plot', '-dpdf', '-r300')


% Create figure with specific size FOR Drift and sigma
figure('Position', [100, 100, 900, 400]);

% Find transition point where kappa becomes >= 1
transition_idx = find(Kappa >= 1, 1, 'first');
if ~isempty(transition_idx)
    transition_eta = eta(transition_idx);
else
    transition_eta = 0.4; % fallback value
end

% Plot data
subplot(1,2,1)
plot(eta, MU, 'b-', 'LineWidth', 2);
hold on
% Add vertical dashed line at transition point
plot([transition_eta transition_eta], ylim, 'k:', 'LineWidth', 1)
hold off

xlabel('\eta', 'FontSize', 14)
ylabel('\mu_{\eta}', 'FontSize', 14)
title('Drift', 'FontSize', 14, 'FontWeight', 'bold')
grid off
xlim([0 1])
ylim([-0.01 0.03])

subplot(1,2,2)
plot(eta, S, 'b-', 'LineWidth', 2);
hold on
% Add horizontal dashed line at y=0
plot([0 1], [0 0], 'k--', 'LineWidth', 1)
% Add vertical dashed line at transition point
plot([transition_eta transition_eta], ylim, 'k:', 'LineWidth', 1)
hold off

xlabel('\eta', 'FontSize', 14)
ylabel('\sigma_{\eta}', 'FontSize', 14)
title('Volatility', 'FontSize', 14, 'FontWeight', 'bold')
grid off
xlim([0 1])
ylim([0 0.1])

% Adjust subplot spacing
set(gcf, 'PaperPositionMode', 'auto')
print('Drift and Vlatility', '-dpdf', '-r300')

