function M = build_M(x, mu, sig)
    N   = length(x);
    dx  = x(2)-x(1);
    dx2 = dx^2;

    dD = -min(mu,0)/dx       + sig.^2/(2*dx2);
    dM = -max(mu,0)/dx + min(mu,0)/dx - sig.^2/dx2;
    dU =  max(mu,0)/dx       + sig.^2/(2*dx2);

    M = spdiags([dD dM dU],[1 0 -1],N,N)';
end