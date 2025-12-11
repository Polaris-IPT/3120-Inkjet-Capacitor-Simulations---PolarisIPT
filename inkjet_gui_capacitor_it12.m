% inkjet_gui_capacitor_zdeflect_fixed.m
% Fixed GUI + 3D animation + voltage staircase for inkjet droplet deflection (z-direction).
% Bug fixes:
%   - moved status text so it doesn't overlap the diameter edit box
%   - added missing L2 parameter and consistent unit handling
%   - fixed cap_end (mm vs m) and time/K calculations
%   - robust start/stop handling and proper in-flight & persistent impacts
%
% Save as inkjet_gui_capacitor_zdeflect_fixed.m and run.

function inkjet_gui_capacitor_it12
    % --------- Default parameters ----------
    params.D_mm  = 3.0;           % paper at x = D_mm (mm)
    params.L1_mm = 0.5;           % capacitor side length and x-length (mm)
    params.cap_start_mm = 1.25;   % cap start x (mm). cap_end = cap_start + L1
    params.L2_mm = 1.25;          % distance from cap exit to paper in mm (kept as param)
    params.W_mm  = 1.0;           % plate separation (mm) -> plates at z = +/- W/2
    params.vx = 10;               % horizontal speed (m/s)
    params.droplet_d = 30e-6;     % droplet diameter (m)
    params.rho = 1000;            % density kg/m^3
    params.q = 1e-12;             % droplet charge (C)
    params.f_fire = 1e4;          % firing frequency Hz
    params.Vmax = 2000;           % max voltage V
    params.dpi = 300;             % resolution
    params.paper_w_mm = 8.5;      % paper width in y (mm)
    params.paper_h_mm = 11.0;     % paper height in z (mm)
    params.animation_speed = 1;   % speed multiplier
    params.markerSize_flight = 80;
    params.markerSize_impact = 60;

    % state
    running = false;
    stopRequested = false;

    % ---------- create UI ----------
    hFig = figure('Name','Inkjet Capacitor GUI (z-deflect) - FIXED','NumberTitle','off',...
        'Units','normalized','Position',[0.02 0.03 0.96 0.9]);

    % left control panel
    ctrlW = 0.22; ctrlH = 0.94;
    hCtrl = uipanel('Parent',hFig,'Title','Controls','FontSize',11,...
        'Units','normalized','Position',[0.01 0.02 ctrlW ctrlH]);

    % 3D axes
    ax3pos = [0.26 0.12 0.5 0.82];
    hAx3 = axes('Parent',hFig,'Units','normalized','Position',ax3pos);
    hold(hAx3,'on'); grid(hAx3,'on'); axis(hAx3,'equal');
    xlabel(hAx3,'x (mm)'); ylabel(hAx3,'y (mm)'); zlabel(hAx3,'z (mm)');
    view(hAx3, [40 20]);
    xlim(hAx3, [0 4]); ylim(hAx3, [-4 4]); zlim(hAx3, [-6 6]);

    % Voltage axes
    axVpos = [0.78 0.12 0.20 0.82];
    hAxV = axes('Parent',hFig,'Units','normalized','Position',axVpos);
    xlabel(hAxV,'time (s)'); ylabel(hAxV,'V (V)');
    title(hAxV,'Voltage staircase V(t)'); grid(hAxV,'on');

    % ---------- Controls layout ----------
    ctrlY = 0.92; step = 0.055;
    addLabel = @(txt,y) uicontrol('Parent',hCtrl,'Style','text','String',txt,...
        'Units','normalized','Position',[0.02 y 0.46 0.04],'HorizontalAlignment','left');
    addEdit = @(val,y,cb) uicontrol('Parent',hCtrl,'Style','edit','String',num2str(val),...
        'Units','normalized','Position',[0.5 y 0.46 0.05],'Callback',cb);

    addLabel('Droplet charge q (C):', ctrlY); hQ = addEdit(params.q, ctrlY, @q_cb);
    ctrlY = ctrlY - step;
    addLabel('Droplet diameter (um):', ctrlY); hDiam = addEdit(params.droplet_d*1e6, ctrlY, @diam_cb);
    ctrlY = ctrlY - step;
    addLabel('Droplet speed vx (m/s):', ctrlY); hVx = addEdit(params.vx, ctrlY, @vx_cb);
    ctrlY = ctrlY - step;
    addLabel('Firing freq f (Hz):', ctrlY); hF = addEdit(params.f_fire, ctrlY, @f_cb);
    ctrlY = ctrlY - step;
    addLabel('Max voltage Vmax (V):', ctrlY); hVmax = addEdit(params.Vmax, ctrlY, @vmax_cb);
    ctrlY = ctrlY - step;
    addLabel('Plate separation W (mm):', ctrlY); hW = addEdit(params.W_mm, ctrlY, @W_cb);
    ctrlY = ctrlY - step;
    addLabel('Capacitor start x (mm):', ctrlY); hCapStart = addEdit(params.cap_start_mm, ctrlY, @capstart_cb);
    ctrlY = ctrlY - step;
    addLabel('L1 (mm):', ctrlY); hL1 = addEdit(params.L1_mm, ctrlY, @L1_cb);
    ctrlY = ctrlY - step;
    addLabel('L2 (mm):', ctrlY); hL2 = addEdit(params.L2_mm, ctrlY, @L2_cb);
    ctrlY = ctrlY - step;
    addLabel('Paper x (mm):', ctrlY); hD = addEdit(params.D_mm, ctrlY, @D_cb);
    ctrlY = ctrlY - step;
    addLabel('Paper width (y) mm:', ctrlY); hPy = addEdit(params.paper_w_mm, ctrlY, @py_cb);
    ctrlY = ctrlY - step;
    addLabel('Paper height (z) mm:', ctrlY); hPz = addEdit(params.paper_h_mm, ctrlY, @pz_cb);
    ctrlY = ctrlY - step;
    addLabel('DPI:', ctrlY); hDPI = addEdit(params.dpi, ctrlY, @dpi_cb);
    ctrlY = ctrlY - step;
    addLabel('Animation speed (mult):', ctrlY); hSpeed = addEdit(params.animation_speed, ctrlY, @speed_cb);
    ctrlY = ctrlY - step;

    % Buttons
    btnH = 0.22;
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Start','Units','normalized',...
        'Position',[0.02 0.05 btnH 0.06],'BackgroundColor',[0.6 1 0.6],'Callback',@start_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Stop','Units','normalized',...
        'Position',[0.27 0.05 btnH 0.06],'BackgroundColor',[1 0.6 0.6],'Callback',@stop_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Reset','Units','normalized',...
        'Position',[0.52 0.05 btnH 0.06],'Callback',@reset_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Compute V(t)','Units','normalized',...
        'Position',[0.02 0.01 0.96 0.035],'Callback',@computeV_cb);

    % Status text moved down to avoid overlap with controls (bottom area above buttons)
    hStatus = uicontrol('Parent',hCtrl,'Style','text','String','Ready','Units','normalized',...
        'Position',[0.02 0.12 0.96 0.045],'BackgroundColor',get(hFig,'Color'),'HorizontalAlignment','left');

    % draw initial scene
    drawScene();

    % ----------------- Functions -----------------
    function drawScene()
        cla(hAx3);
        hold(hAx3,'on'); grid(hAx3,'on'); axis(hAx3,'equal');
        xlabel(hAx3,'x (mm)'); ylabel(hAx3,'y (mm)'); zlabel(hAx3,'z (mm)');
        xlim(hAx3, [0 4]); ylim(hAx3, [-4 4]); zlim(hAx3, [-6 6]);
        view(hAx3,[40 20]);

        % nozzle point
        scatter3(hAx3, 0, 0, 0, 140, 'k', 'filled');
        text(hAx3, 0.12, -0.25, 0, 'Nozzle', 'FontSize', 10);

        % capacitor plates (square in xy-plane at z = +/- W/2)
        x_cap_start_mm = params.cap_start_mm;
        x_cap_end_mm = x_cap_start_mm + params.L1_mm;
        y_plate_min = -params.L1_mm/2;
        y_plate_max =  params.L1_mm/2;
        z_top_mm = params.W_mm/2;
        z_bot_mm = -params.W_mm/2;

        px = [x_cap_start_mm, x_cap_end_mm, x_cap_end_mm, x_cap_start_mm];
        py = [y_plate_min, y_plate_min, y_plate_max, y_plate_max];
        pz_top = z_top_mm * ones(1,4);
        pz_bot = z_bot_mm * ones(1,4);

        patch(hAx3, px, py, pz_top, [0.8 0.9 1], 'FaceAlpha',0.65,'EdgeColor','none');
        patch(hAx3, px, py, pz_bot, [0.8 0.9 1], 'FaceAlpha',0.65,'EdgeColor','none');

        % plate edges
        plot3(hAx3, [x_cap_start_mm x_cap_end_mm],[y_plate_min y_plate_min],[z_top_mm z_top_mm],'k-','LineWidth',0.8);
        plot3(hAx3, [x_cap_start_mm x_cap_end_mm],[y_plate_max y_plate_max],[z_top_mm z_top_mm],'k-','LineWidth',0.8);
        plot3(hAx3, [x_cap_start_mm x_cap_end_mm],[y_plate_min y_plate_min],[z_bot_mm z_bot_mm],'k-','LineWidth',0.8);
        plot3(hAx3, [x_cap_start_mm x_cap_end_mm],[y_plate_max y_plate_max],[z_bot_mm z_bot_mm],'k-','LineWidth',0.8);

        % paper at x = D_mm (yz plane)
        Dmm = params.D_mm;
        paper_y = [ params.paper_w_mm/2, -params.paper_w_mm/2, -params.paper_w_mm/2, params.paper_w_mm/2 ];
        paper_z = [ params.paper_h_mm/2,  params.paper_h_mm/2, -params.paper_h_mm/2, -params.paper_h_mm/2 ];
        px_vals = Dmm * ones(1,4);
        patch(hAx3, px_vals, paper_y, paper_z, [1 1 0.9], 'FaceAlpha',0.95,'EdgeColor','k');

        scatter3(hAx3, Dmm, 0, 0, 80, 'filled', 'MarkerEdgeColor','k');
        title(hAx3, '3D View: Nozzle -> Capacitor Plates -> Paper');
        drawnow;
    end

    function computeV_cb(~,~)
        compute_and_plot_voltage();
    end

    function compute_and_plot_voltage()
        mm2m = 1e-3;
        Dm = params.D_mm * mm2m;
        L1m = params.L1_mm * mm2m;
        Wm = params.W_mm * mm2m;
        cap_start_m = params.cap_start_mm * mm2m;
        cap_end_m = (params.cap_start_mm + params.L1_mm) * mm2m;

        r = params.droplet_d / 2;
        m = (4/3) * pi * r^3 * params.rho;

        T_inside = L1m / params.vx;
        t_after = (Dm - cap_end_m) / params.vx;

        dot_spacing = (25.4e-3) / params.dpi;
        paper_h_m = params.paper_h_mm * mm2m;
        N_dots = max(1, round(paper_h_m / dot_spacing));
        z_targets = linspace(paper_h_m/2, -paper_h_m/2, N_dots);

        K = (params.q / (m * Wm)) * (0.5 * T_inside^2 + T_inside * t_after);
        V_required = z_targets / K;
        V_clip = V_required;
        V_clip(abs(V_clip) > params.Vmax) = sign(V_clip(abs(V_clip) > params.Vmax)) .* params.Vmax;

        % times
        t_fire = (0:N_dots-1) / params.f_fire;
        t_exit = t_fire + cap_end_m / params.vx;
        t_plot = [t_exit, t_exit(end) + 1/params.f_fire];
        V_plot = [V_clip, V_clip(end)];

        axes(hAxV);
        cla(hAxV);
        stairs(t_plot, V_plot, 'LineWidth', 1.5);
        xlabel(hAxV,'time (s)'); ylabel(hAxV,'V (V)');
        title(hAxV,'Voltage staircase V(t)');
        grid(hAxV,'on');
    end

    function start_cb(~,~)
        if running
            return;
        end
        running = true;
        stopRequested = false;
        set(hStatus,'String','Running');
        drawnow;
        compute_and_run_simulation();
    end

    function stop_cb(~,~)
        stopRequested = true;
        running = false;
        set(hStatus,'String','Stop requested');
        drawnow;
    end

    function reset_cb(~,~)
        stopRequested = true;
        running = false;
        set(hStatus,'String','Reset - ready');
        drawScene();
        axes(hAxV); cla(hAxV);
        drawnow;
    end

    % ----------- param edit callbacks (update params only) ----------
    function q_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.q = val; end; end
    function diam_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.droplet_d = val*1e-6; end; end
    function vx_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.vx = val; end; end
    function f_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.f_fire = val; end; end
    function vmax_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.Vmax = val; end; end
    function W_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.W_mm = val; end; end
    function capstart_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.cap_start_mm = val; end; end
    function L1_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.L1_mm = val; end; end
    function L2_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.L2_mm = val; end; end
    function D_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.D_mm = val; end; end
    function py_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.paper_w_mm = val; end; end
    function pz_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.paper_h_mm = val; end; end
    function dpi_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.dpi = val; end; end
    function speed_cb(hObj,~); val = str2double(get(hObj,'String')); if ~isnan(val); params.animation_speed = max(0.01,val); end; end

    % ----------- core simulation (fixed) ----------------
    function compute_and_run_simulation()
        mm2m = 1e-3;
        Dm = params.D_mm * mm2m;
        L1m = params.L1_mm * mm2m;
        Wm = params.W_mm * mm2m;
        cap_start_m = params.cap_start_mm * mm2m;
        cap_end_m = (params.cap_start_mm + params.L1_mm) * mm2m;

        r = params.droplet_d / 2;
        m = (4/3) * pi * r^3 * params.rho;

        T_inside = L1m / params.vx;
        t_after = (Dm - cap_end_m) / params.vx;  % time after cap exit to paper
        t_total = Dm / params.vx;

        % dot spacing & targets
        dot_spacing = (25.4e-3) / params.dpi;
        paper_h_m = params.paper_h_mm * mm2m;
        N_dots = max(1, round(paper_h_m / dot_spacing));
        z_targets = linspace(paper_h_m/2, -paper_h_m/2, N_dots);

        % K factor and voltages
        K = (params.q / (m * Wm)) * (0.5 * T_inside^2 + T_inside * t_after);
        V_required = z_targets / K;
        V_required(abs(V_required) > params.Vmax) = sign(V_required(abs(V_required) > params.Vmax)) .* params.Vmax;

        % times for staircase
        t_fire = (0:N_dots-1) / params.f_fire;
        t_exit = t_fire + cap_end_m / params.vx; % seconds
        t_impact = t_fire + t_total;

        % plot staircase
        axes(hAxV); cla(hAxV);
        t_plot = [t_exit, t_exit(end) + 1/params.f_fire];
        V_plot = [V_required, V_required(end)];
        stairs(hAxV, t_plot, V_plot, 'LineWidth', 1.5);
        xlabel(hAxV,'time (s)'); ylabel(hAxV,'V (V)'); title(hAxV,'Voltage staircase V(t)');
        grid(hAxV,'on');

        % prepare scene
        drawScene();
        cmap = feval('jet', N_dots);

        % persistent impact markers container
        impactMarkersLocal = gobjects(N_dots,1);

        % in-flight marker handle
        hFlight = scatter3(hAx3, nan, nan, nan, params.markerSize_flight, 'filled','MarkerEdgeColor','k');

        % convert cap positions to meters for physics
        x_cap_start_m = cap_start_m;
        x_cap_end_m = cap_end_m;

        set(hStatus,'String','Simulating...'); drawnow;

        % Sequential droplets (visual clarity). If you prefer simultaneous,
        % I can switch to a time-stepping multi-droplet mode.
        for k = 1:N_dots
            if stopRequested
                stopRequested = false;
                break;
            end

            Vk = V_required(k);
            Ez = Vk / Wm;                      % V/m along z
            az_inside = (params.q * Ez) / m;   % m/s^2 while inside cap

            % init conditions SI
            x = 0; z = 0; y = 0;
            vx_now = params.vx; vz_now = 0;

            % integration steps chosen for visual smoothness
            steps = max(40, round(100 / params.animation_speed));
            dt_local = t_total / steps;

            for s = 1:steps
                % determine if inside capacitor region (use current x)
                if x >= x_cap_start_m && x <= x_cap_end_m
                    az = az_inside;
                else
                    az = 0;
                end

                % integrate
                vz_now = vz_now + az * dt_local;
                x = x + vx_now * dt_local;
                z = z + vz_now * dt_local;

                % update flight marker (in mm)
                set(hFlight, 'XData', x/mm2m, 'YData', y/mm2m, 'ZData', z/mm2m, ...
                    'CData', cmap(k,:), 'SizeData', params.markerSize_flight);
                drawnow limitrate;
            end

            % place persistent impact marker at paper x = params.D_mm
            landing_y_mm = y / mm2m;
            landing_z_mm = z / mm2m;
            impactMarkersLocal(k) = scatter3(hAx3, params.D_mm, landing_y_mm, landing_z_mm, ...
                params.markerSize_impact, 'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor','k');

            % hide in-flight marker for next cycle
            set(hFlight,'XData',nan,'YData',nan,'ZData',nan);
            pause(max(0, 0.002 / params.animation_speed));
        end

        set(hStatus,'String','Simulation complete');
        running = false;
    end
end
