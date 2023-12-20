function movies()
    % change for different experiments 1, 2, or 3
    mov_num = 1;
    % change for different plots 1 or 2
    plt_mode = 1;
    level = 8;
    % set to 1 to make movie
    avienable = 0;

    % Initial Conditions
    if mov_num == 1
        tmax = 0.04;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.3, 0.5, 0.055, 0.055, 20, 10];
        vtype = 1;
        vpar = [0.5, 0.7, 0.2, 0.8, 10^6];
        plot_title = ["Scattering Off a Rectangular Barrier", ...
            "(x_{min}=0.5, x_{max}=0.7, y_{min}=0.2, y_{max}=0.8, V_{c}=10^{6})"];
    elseif mov_num == 2
        tmax = 0.04;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.3, 0.5, 0.055, 0.055, 0, 0];
        vtype = 1;
        vpar = [0.5, 0.7, 0.2, 0.8, -10^4];
        plot_title = ["Scattering Off a Rectangular Well", ...
            "(x_{min}=0.5, x_{max}=0.7, y_{min}=0.2, y_{max}=0.8, V_{c}=-10^{4})"];
    elseif mov_num == 3
        tmax = 0.04;
        lambda = 0.01;
        idtype = 1;
        idpar = [0.5, 0.1, 0.055, 0.055, 0, 15];
        vtype = 2;
        vpar = [0.3, 0.4, 0.6, 0.7, 10^6];
        plot_title = [" Scattering Through a Double Slit", ...
            "(x_{1}=0.3, x_{2}=0.4, x_{3}=0.6, x_{4}=0.7, V_{c}=10^{6})"];
    end

    % calculate solution
    [x, y, t, psi, psire, psiim, psimod, v] = ...
        sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);
    
    aspect_ratio = [1 1 1];
    % for dowesampling solution
    ds = 4;
    avifilename = 'test_movie.avi';

    % axis
    x_lim = [0 1];
    y_lim = [0 1];
    z_lim = [0 1];

    % get avi object
    if avienable
        aviobj = VideoWriter(avifilename);
        open(aviobj);
    end

    for ii = 1:size(t, 2)
        if plt_mode == 1
            % downsample
            psimod_ds = psimod(ii, 1:ds:end, 1:ds:end);
            % made 3d surface plot
            surf(x(1:ds:end), y(1:ds:end), squeeze(psimod_ds));
            shading faceted;
            view(-37.5, 30);
            %set axis
            xlim(x_lim);
            ylim(y_lim);
            zlim(z_lim);
            % set aspect ratio
            pbaspect(aspect_ratio);
            %set labels
            xlabel("x", "interpreter", "tex")
            ylabel("y", "interpreter", "tex")
            zlabel("|\psi(x, y, t)|", "interpreter", "tex")
        elseif plt_mode == 2
            % make contour plot
            contourf(x, y, squeeze(psimod(ii, :, :)));
            % set labels
            xlabel("x", "interpreter", "tex")
            ylabel("y", "interpreter", "tex")
        end

        title(plot_title, "interpreter", "tex");
        drawnow;
        if avienable
            if t == 0
                framecount = aviframerate*5;
            else
                framecount = 1;
            end
            for iframe = 1:framecount
                writeVideo(aviobj, getframe(gcf));
            end
        end

        if mod(ii, 10) == 0
            fprintf("Iteration " + string(ii) + " out of " + string(size(t, 2)) + "\n");
        end
    end
    
    if avienable
        close(aviobj);
        fprintf('Created video file: %s\n', avifilename);
    end

 end 