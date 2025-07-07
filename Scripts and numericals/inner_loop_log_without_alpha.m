function [Q, SSQ, Kappa, Iota] = inner_loop_log_without_alpha(eta, q0, a_e, a_h, rho_e, rho_h, sigma, phi)
    N    = length(eta);
    deta = [eta(1); diff(eta)];  % imposes the correct grid step for numerical derivative at η^e = 0

    % variables
    Q     = ones(N,1);    % price of capital q
    SSQ   = zeros(N,1);   % σ + σ^q
    Kappa = zeros(N,1);   % capital fraction of experts κ

    Rho   = eta * rho_e + (1 - eta) * rho_h;  % average consumption-to-networth ratio

    % Initiate the loop
    kappa  = 0;
    q_old  = q0;
    q      = q0;
    ssq    = sigma;

    % Iterate over eta
    % At each step apply Newton’s method to F(z) = 0 where z = [q, kappa, ssq]’
    % Use chi = kappa
    for i = 1:N
        % Compute F(z_{n-1})
        F = [
            kappa * (a_e - a_h) + a_h - (q - 1)/phi - q * Rho(i);
            ssq * (q - (q - q_old)/deta(i) * (kappa - eta(i))) - sigma * q;
            a_e - a_h - q *(kappa - eta(i)) / (eta(i)*(1 - eta(i))) * ssq^2
        ];

        % Construct Jacobian J^{n-1}
        J = zeros(3,3);
        J(1,:) = [-1/phi - Rho(i), a_e - a_h, 0];
        J(2,:) = [
            ssq * (1 - (kappa - eta(i))/deta(i)) - sigma, ...
            -ssq * (q - q_old)/deta(i), ...
            q - (q - q_old)/deta(i) * (kappa - eta(i))
        ];
        J(3,:) = [
            -(kappa - eta(i)) / (eta(i)*(1 - eta(i))) * ssq^2, ...
            -q/ (eta(i)*(1 - eta(i))) * ssq^2, ...
            -2 * q * (kappa - eta(i)) / (eta(i)*(1 - eta(i))) * ssq
        ];

        % Iterate, obtain z_{n}
        z = [q; kappa; ssq] - J \ F;

        % If the new kappa is larger than 1, break
        if z(2) >= 1
            break;
        end

        % Update variables
        q     = z(1);
        kappa = z(2);
        ssq   = z(3);

        % save results
        Q(i)     = q;
        Kappa(i) = kappa;
        SSQ(i)   = ssq;
        q_old    = q;
    end

    % Set kappa = 1, use chi = kappa and compute the rest
    n1 = i;
    for i = n1:N
        q  = (1 + a_e * phi) / (1 + Rho(i) * phi);
        qp = (q - q_old) / deta(i);

        Q(i)      = q;
        Kappa(i)  = 1;
        SSQ(i)    = sigma / (1 - (1 - eta(i)) * qp / q);
        q_old     = q;
    end

    Iota = (Q - 1) / phi;
end
