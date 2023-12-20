function ctest_2d()
    % get solutions for convtest for different levels
    mx = 2;
    my = 3;
    [x6 y6 t6 psi6 psire6 psiim6 psimod6 v6] = ...
        sch_2d_adi(0.05, 6, 0.05, 0, [mx, my], 0, []);
    [x7 y7 t7 psi7 psire7 psiim7 psimod7 v7] = ...
        sch_2d_adi(0.05, 7, 0.05, 0, [mx, my], 0, []);
    [x8 y8 t8 psi8 psire8 psiim8 psimod8 v8] = ...
        sch_2d_adi(0.05, 8, 0.05, 0, [mx, my], 0, []);
    [x9 y9 t9 psi9 psire9 psiim9 psimod9 v9] = ...
        sch_2d_adi(0.05, 9, 0.05, 0, [mx, my], 0, []);

    % test rectangular barrier potential works
    sch_2d_adi(0.25, 6, 0.1, 1, [2, 3, 4, 5, 6, 7], 1, [1 2 3 4 5]);

    % test double slit potential works
    sch_2d_adi(0.25, 6, 0.1, 1, [2, 3, 4, 5, 6, 7], 2, [1 2 3 4 5]);

    fprintf('Done Running sch_2d_adi for each case, now plotting \n');
    
    % coarsen solutions to match size of lmin solution
    psi7 = psi7(1:2:end, 1:2:end, 1:2:end);
    psi8 = psi8(1:4:end, 1:4:end, 1:4:end);
    psi9 = psi9(1:8:end, 1:8:end, 1:8:end);

    dpsi6 = psi7 - psi6;
    dpsi7 = psi8 - psi7;
    dpsi8 = psi9 - psi8;

    % calculate rms values of dpsi
    norm_dpsi6 = psi_norm(dpsi6);
    norm_dpsi7 = psi_norm(dpsi7);
    norm_dpsi8 = psi_norm(dpsi8);

    close all;
    figure;
    hold on;
    titlestr = sprintf('Scaled Spatial Norm vs Time');
    title(titlestr, 'interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', [0.25, 0.42, 0.31]);
    xlabel('Time');
    ylabel('Spatial Norm');
    plot(t6, norm_dpsi6, 'r-.o');
    plot(t6, 4 * norm_dpsi7, 'g-.+'); 
    plot(t6, 4^2 * norm_dpsi8, 'b-.*');
    legend('||dpsi^{6}||_2', '4*||dpsi^{7}||_2', '16*||dpsi^{8}||_2');
    hold off;

    % calculate exact solution
    psi_exact = zeros(size(t6, 2), size(x6, 2), size(y6, 2));
    for idt = 1:size(t6, 2)
        for idx = 1:size(x6, 2)
            psi_exact(idt, idx, :) = exp(-1i*(mx^2+my^2)*pi^2*t6(idt)) * ...
                sin(mx*pi*x6(idx)) * sin(my*pi*y6);
        end
    end

    % calculate rms values of E
    norm_E6 = psi_norm(psi_exact - psi6);
    norm_E7 = psi_norm(psi_exact - psi7);
    norm_E8 = psi_norm(psi_exact - psi8);
    norm_E9 = psi_norm(psi_exact - psi9);

    figure;
    hold on;
    titlestr = sprintf('Scaled Spatial Norm vs Time');
    title(titlestr, 'interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', [0.25, 0.42, 0.31]);
    xlabel('Time');
    ylabel('Spatial Norm');
    plot(t6, norm_E6, 'r-.o');
    plot(t6, 4 * norm_E7, 'g-.+'); 
    plot(t6, 4^2 * norm_E8, 'b-.*');
    plot(t6, 4^3 * norm_E9, 'c-.*');
    legend('||E(psi^{6})||_2', '4*||E(psi^{7})||_2', '16*||E(psi^{8})||_2', '64*||E(psi^{9})||_2');
    hold off;

end

% Output
% norm: Vector of norms [nt]
function norm = psi_norm(psi)
    [nt, nx, ny] = size(psi);
    norm = zeros(nt, 1);
    for n = 1:nt
        for ii = 1:nx
            for jj = 1:ny
                norm(n) = norm(n) + abs(psi(n, ii, jj))^2;
            end
        end
    end
    norm = sqrt(norm / nx / ny);
end 