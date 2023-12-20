% Inputs
%
% tmax: Maximum integration time
% level: Discretization level
% lambda: dt/dx
% idtype: Selects initial data type
% idpar: Vector of initial data parameters
% vtype: Selects potential type
% vpar: Vector of potential parameters
%
% Outputs
%
% x: Vector of x coordinates [nx]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx]
% psire Array of computed psi_re values [nt x nx]
% psiim Array of computed psi_im values [nt x nx]
% psimod Array of computed sqrt(psi psi*) values [nt x nx]
% prob Array of computed running integral values [nt x nx]
% v Array of potential values [nx]
function [x t psi psire psiim psimod prob v] = ...
    sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar)

    % Discretization Parameters
    nx = 2^level + 1;
    dx = 2^(-level);
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    x = linspace(0.0, 1.0, nx);
    t = linspace(0.0, tmax, nt);

    % Solution Parameters
    psi = zeros(nt, nx);
    psire = zeros(nt, nx);
    psiim = zeros(nt, nx);
    psimod = zeros(nt, nx);
    prob = zeros(nt, nx);

    % Initial Data
    if idtype == 0
        psi(1, :) = sin(idpar(1) * pi * x);
    elseif idtype == 1
        psi(1, :) = exp(1i*idpar(3)*x) .* exp(-((x-idpar(1)) / idpar(2)).^2);
    else
        fprintf("Invalid idtype");
        return
    end
    
    % Potential Data
    if vtype == 0
        v = zeros(nx, 1);
    elseif vtype == 1
        v = zeros(nx, 1);
        for idx = 1:nx
            if vpar(1) <= x(idx) && x(idx) <= vpar(2)
                v(idx) = vpar(3);
            end
        end
    else
        fprintf("Invalid idtype");
        return
    end
    
    % Initialize storage for sparse matrix and RHS
    dl = zeros(nx,1);
    d  = zeros(nx,1);
    du = zeros(nx,1);
    f  = zeros(nx,1);

    % Set up tridiagonal system
    dl = 0.5 / dx^2 .* ones(nx, 1);
    d  = (1i / dt - 1.0 / dx^2 - v / 2) .* ones(nx,1);
    du = dl;

    % Fix up boundary cases
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx-1) = 0.0;
    d(nx) = 1.0;

    % Define sparse matrix
    A = spdiags([dl d du], -1:1, nx, nx);

    % Compute solution using CN scheme
    for n = 1 : nt-1
        % Define RHS of linear system
        f(2:nx-1) = 1i * psi(n, 2:nx-1) / dt - 0.5 * ( ...
            psi(n, 1:nx-2) - 2 * psi(n, 2:nx-1) + psi(n, 3:nx)) / dx^2 ...
            + 0.5 * v(2:nx-1).' .* psi(n, 2:nx-1);
        f(1) = 0.0;
        f(nx) = 0.0;
        % Solve system, thus updating approximation to next time step
        psi(n+1, :) = A \ f;
    end

    psire = real(psi);
    psiim = imag(psi);
    psimod = sqrt(psi .* conj(psi));
    prob = zeros(nt, nx);
    
    % calculate running integral of the probability density using trapezoidal approximation
    ro = psimod.^2;
    for ii = 1 : nt
        sum = 0;
        for jj = 2 : nx
            sum = sum + 0.5 * (ro(ii, jj-1) + ro(ii,jj)) * (x(jj) - x(jj-1));
            prob(ii, jj) = sum;
        end
    end
end