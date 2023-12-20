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
% y: Vector of y coordinates [ny]
% t: Vector of t coordinates [nt]
% psi: Array of computed psi values [nt x nx x ny]
% psire Array of computed psi_re values [nt x nx x ny]
% psiim Array of computed psi_im values [nt x nx x ny]
% psimod Array of computed sqrt(psi psi*) values [nt x nx x ny]
% v Array of potential values [nx x ny]
function [x y t psi psire psiim psimod v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar)

    % Discretization Parameters
    nx = 2^level + 1;
    ny = nx;
    dx = 2^(-level);
    dy = dx;
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;
    x = linspace(0.0, 1.0, nx);
    y = linspace(0.0, 1.0, ny);
    t = linspace(0.0, tmax, nt);

    % Solution Parameters
    psi = zeros(nt, nx, ny);
    psire = zeros(nt, nx, ny);
    psiim = zeros(nt, nx, ny);
    psimod = zeros(nt, nx, ny);

    % Initial Data
    if idtype == 0
        for idx = 1:nx
            psi(1, idx, :) = sin(idpar(1) * pi * x(idx)) * sin(idpar(2) * pi * y);
        end
    elseif idtype == 1
        for idx = 1:nx
            psi(1, idx, :) = exp(1i*idpar(5)*x(idx)) * ...
                exp(1i*idpar(6)*y) .* ...
                exp(-((x(idx) - idpar(1))^2 / idpar(3)^2 +...
                (y - idpar(2)).^2 / idpar(4)^2));
        end
    else
        fprintf("Invalid idtype");
        return
    end
    
    % Potential Data
    if vtype == 0
        v = zeros(nx, ny);
    elseif vtype == 1
        v = zeros(nx, ny);
        for ii = 1:nx
            for jj = 1:ny
                if vpar(1) <= x(ii) && x(ii) <= vpar(2) && ...
                        vpar(3) <= y(jj) && y(jj) <= vpar(4)
                    v(ii, jj) = vpar(5);
                end
            end
        end
    elseif vtype == 2
        v = zeros(nx, ny);
        jp = (ny-1)/4 + 1;
        v(:, jp) = vpar(5);
        v(:, jp+1) = vpar(5);
        for idx = 1:nx
            if (vpar(1) <= x(idx) && x(idx) <= vpar(2)) || ...
                    (vpar(3) <= x(idx) && x(idx) <= vpar(4))
                v(idx, jp) = 0;
                v(idx, jp+1) = 0;
            end
        end
    else
        fprintf("Invalid idtype");
        return
    end
    
    % Initialize storage for sparse matrix 1,2 and RHS 1,2
    dl = zeros(nx,1);
    d  = zeros(nx,1);
    du = zeros(nx,1);
    f1  = zeros(nx,1);
    f2  = zeros(ny,1);

    % Set up tridiagonal system 1
    dl = -0.5i * dt / dx^2 * ones(nx, 1);
    d  = (1 + 1i * dt / dx^2) .* ones(nx,1);
    du = dl;

    % Fix up boundary cases
    d(1) = 1.0;
    du(2) = 0.0;
    dl(nx-1) = 0.0;
    d(nx) = 1.0;

    % Define sparse matrix 1
    A_half = spdiags([dl d du], -1:1, nx, nx);

    % Set up tridiagonal system 2 in the loop
    % dl and du unchanged because ny=nx
    % d is x and y dept. so add potential in the loop

    % Compute solution using ADI scheme
    for n = 1 : nt-1
        % Solve half step system
        psi_half = zeros(nx, ny);
        for k = 2:ny-1
            % Define RHS of linear system
            dyy = (psi(n,:,k+1) - 2*psi(n,:,k) + psi(n,:,k-1)) / dy^2;
            % right bracket calculation (intermediate value)
            psi_inter = psi(n,:,k) + 0.5i*dt*dyy - 0.5i*dt*v(:,k).'.*psi(n,:,k);

            % full RHS calculation
            dxx = (psi_inter(1:nx-2) - 2*psi_inter(2:nx-1) + psi_inter(3:nx)) / dx^2;
            f1(2:nx-1) = psi_inter(2:nx-1) + 0.5i*dt*dxx;

            % BCs
            f1(1) = 0.0;
            f1(ny) = 0.0;

            % Solve system 1, thus updating approximation to next time step
            psi_half(:, k) = A_half \ f1;
        end

        % Solve one step system
        for h = 2:nx-1

            % add potential to sparse matrix
            dv = d + 0.5i*dt*v(h,:).';
            % Fix up boundary cases
            dv(1) = 1.0;
            dv(nx) = 1.0;
            % Define sparse matrix 2 with potential
            A_one = spdiags([dl dv du], -1:1, ny, ny);

            % Define RHS of linear system
            f2 = psi_half(h, :).';

            % BCs
            f2(1) = 0.0;
            f2(nx) = 0.0;

            % Solve system 2, thus updating approximation to next time step
            psi(n+1, h, :) = A_one \ f2;
        end
    end

    psire = real(psi);
    psiim = imag(psi);
    psimod = sqrt(psi .* conj(psi));
end