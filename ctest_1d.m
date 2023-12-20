function ctest_1d()
    % get solutions for convtest 1 for different levels
    m = 3;
    [x1_6 t1_6 psi1_6 psire1_6 psiim1_6 psimod1_6 prob1_6 v1_6] = ...
        sch_1d_cn(0.25, 6, 0.1, 0, [m], 0, []);
    [x1_7 t1_7 psi1_7 psire1_7 psiim1_7 psimod1_7 prob1_7 v1_7] = ...
        sch_1d_cn(0.25, 7, 0.1, 0, [m], 0, []);
    [x1_8 t1_8 psi1_8 psire1_8 psiim1_8 psimod1_8 prob1_8 v1_8] = ...
        sch_1d_cn(0.25, 8, 0.1, 0, [m], 0, []);
    [x1_9 t1_9 psi1_9 psire1_9 psiim1_9 psimod1_9 prob1_9 v1_9] = ...
        sch_1d_cn(0.25, 9, 0.1, 0, [m], 0, []);
    %[x2 t2 psi2 psire2 psiim2 psimod2 prob2 v2] = sch_1d_cn(0.01, 4, 0.01, 1, [0.50 0.075 0.0], 0, []);

    % coarsen solutions to match size of lmin solution
    psi1_7 = psi1_7(1:2:end, 1:2:end);
    psi1_8 = psi1_8(1:4:end, 1:4:end);
    psi1_9 = psi1_9(1:8:end, 1:8:end);

    dpsi6 = psi1_7 - psi1_6;
    dpsi7 = psi1_8 - psi1_7;
    dpsi8 = psi1_9 - psi1_8;

    % calculate rms values of dpsi
    rms_dpsi6 = rms(dpsi6, 2);
    rms_dpsi7 = rms(dpsi7, 2);
    rms_dpsi8 = rms(dpsi8, 2);

    close all;
    figure;
    hold on;
    titlestr = sprintf('Scaled Spatial Norm vs Time');
    title(titlestr, 'interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', [0.25, 0.42, 0.31]);
    xlabel('Time');
    ylabel('Spatial Norm');
    plot(t1_6, rms_dpsi6, 'r-.o');
    plot(t1_6, 4 * rms_dpsi7, 'g-.+'); 
    plot(t1_6, 4^2 * rms_dpsi8, 'b-.*');
    legend('||dpsi^{6}||_2', '4*||dpsi^{7}||_2', '16*||dpsi^{8}||_2');
    hold off;

    % calculate exact solution
    psi_exact = zeros(size(t1_6, 2), size(x1_6, 2));
    for idx = 1:size(t1_6, 2)
        psi_exact(idx, :) = exp(-1i * m^2 * pi^2 * t1_6(idx)) * sin(m * pi * x1_6);
    end

    % calculate rms values of E
    rms_E6 = rms(psi_exact - psi1_6, 2);
    rms_E7 = rms(psi_exact - psi1_7, 2);
    rms_E8 = rms(psi_exact - psi1_8, 2);
    rms_E9 = rms(psi_exact - psi1_9, 2);

    figure;
    hold on;
    titlestr = sprintf('Scaled Spatial Norm vs Time');
    title(titlestr, 'interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', [0.25, 0.42, 0.31]);
    xlabel('Time');
    ylabel('Spatial Norm');
    plot(t1_6, rms_E6, 'r-.o');
    plot(t1_6, 4 * rms_E7, 'g-.+'); 
    plot(t1_6, 4^2 * rms_E8, 'b-.*');
    plot(t1_6, 4^3 * rms_E9, 'c-.*');
    legend('||E(psi^{6})||_2', '4*||E(psi^{7})||_2', '16*||E(psi^{8})||_2', '64*||E(psi^{9})||_2');
    hold off;
    
end