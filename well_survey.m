function well_survey()
    % experiment variables
    level = 9;
    tmax = 0.10;
    lambda = 0.01;
    x1 = 0.6;
    x2 = 0.8;
    num_exp = 251;
    lnV0 = linspace(2, 10, num_exp);
    V0 = exp(lnV0);

    % size variables
    nx = 2^level + 1;
    dx = 2^(-level);
    dt = lambda * dx;
    nt = round(tmax / dt) + 1;

    % calculate output vectors
    Fe = zeros(num_exp, 1);

    % run experiments
    for idx = 1:num_exp
        [x t psi psire psiim psimod prob v] = ...
            sch_1d_cn(tmax, level, lambda, 1, [0.40, 0.075, 0.0], 1, [0.6 0.8 -V0(idx)]);

        % p_j^n = prob[x_j, t^n] is the convention
        % calculate normalized temporal average
        temp_avg_prob = mean(prob);
        temp_avg_prob = temp_avg_prob / temp_avg_prob(nx);

        % find closest points in x to x1 and x2
        [~, x1_idx] = min(abs(x - x1));
        [~, x2_idx] = min(abs(x - x2));

        % calculate fractional probability
        Fe(idx) = (temp_avg_prob(x2_idx) - temp_avg_prob(x1_idx)) / (x2 - x1);
        if Fe >= 1
            fprintf("Broken Rule");
        end

        if mod(idx, 30) == 0
            fprintf("Iteration " + string(idx) + "\n");
        end
    end

    close all;
    figure;
    hold on;
    titlestr = sprintf('Well Survey');
    title(titlestr, 'interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold', ...
        'Color', [0.25, 0.42, 0.31]);
    xlabel('ln(|V0|))');
    ylabel('ln(F_{e}(' + string(x1) + ', ' + string(x2) + '))');
    plot(lnV0, log(Fe), 'r-.o');
    hold off;
end