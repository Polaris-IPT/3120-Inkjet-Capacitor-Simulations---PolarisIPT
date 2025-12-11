function inkjet_gui_capacitor_zdeflect_fixed_timed
% inkjet_gui_capacitor_zdeflect_fixed_timed.m
% Vertical-only GUI (I-shape) with timing shown in GUI.
% Save and run this file in MATLAB.

    % ---------------- Default parameters ----------------
    params.D_mm  = 3.0;
    params.L1_mm = 0.5;
    params.cap_start_mm = 1.25;
    params.L2_mm = 1.25;
    params.W_mm  = 1.0;
    params.vx = 10;                 % m/s
    params.droplet_d = 30e-6;       % m
    params.rho = 1000;              % kg/m^3
    params.q = 1e-12;               % C
    params.f_fire = 1e4;            % Hz
    params.Vmax = 2000;             % V
    params.dpi = 300;
    params.paper_w_mm = 8.5;
    params.paper_h_mm = 11.0;
    params.animation_speed = 1;
    params.markerSize_flight = 80;
    params.markerSize_impact = 60;

    running = false;
    stopRequested = false;

    % ---------------- Create GUI ----------------
    hFig = figure('Name','Inkjet Vertical GUI (I) - timed','NumberTitle','off',...
        'Units','normalized','Position',[0.02 0.03 0.96 0.9]);

    % Left control panel
    ctrlW = 0.22; ctrlH = 0.94;
    hCtrl = uipanel('Parent',hFig,'Title','Controls','FontSize',11,...
        'Units','normalized','Position',[0.01 0.02 ctrlW ctrlH]);

    % 3D axes
    ax3pos = [0.26 0.12 0.5 0.82];
    hAx3 = axes('Parent',hFig,'Units','normalized','Position',ax3pos);
    hold(hAx3,'on'); grid(hAx3,'on'); axis(hAx3,'equal');
    xlabel(hAx3,'x (mm)'); ylabel(hAx3,'y (mm)'); zlabel(hAx3,'z (mm)');
    view(hAx3,[40 20]);
    xlim(hAx3,[0 4]); ylim(hAx3,[-4 4]); zlim(hAx3,[-6 6]);

    % Voltage plot axis (single)
    axVpos = [0.78 0.12 0.20 0.82];
    hAxV = axes('Parent',hFig,'Units','normalized','Position',axVpos);
    xlabel(hAxV,'time (s)'); ylabel(hAxV,'V (V)');
    title(hAxV,'Voltage staircase V_z(t)'); grid(hAxV,'on');

    % Controls layout
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

    % Buttons & status
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Start','Units','normalized',...
        'Position',[0.02 0.05 0.22 0.06],'BackgroundColor',[0.6 1 0.6],'Callback',@start_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Stop','Units','normalized',...
        'Position',[0.27 0.05 0.22 0.06],'BackgroundColor',[1 0.6 0.6],'Callback',@stop_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Reset','Units','normalized',...
        'Position',[0.52 0.05 0.22 0.06],'Callback',@reset_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Compute V(t)','Units','normalized',...
        'Position',[0.02 0.01 0.96 0.035],'Callback',@computeV_cb);

    % Status and timing displays (GUI)
    hStatus = uicontrol('Parent',hCtrl,'Style','text','String','Ready','Units','normalized',...
        'Position',[0.02 0.12 0.96 0.045],'BackgroundColor',get(hFig,'Color'),'HorizontalAlignment','left');

    % Timing displays
    addLabel('First droplet time (s):', 0.18);
    hFirstTime = uicontrol('Parent',hCtrl,'Style','text','String','N/A','Units','normalized',...
        'Position',[0.5 0.18 0.46 0.04],'HorizontalAlignment','left');
    addLabel('Total I draw time (s):', 0.13);
    hTotalTime = uicontrol('Parent',hCtrl,'Style','text','String','N/A','Units','normalized',...
        'Position',[0.5 0.13 0.46 0.04],'HorizontalAlignment','left');

    % Initial scene
    drawScene();

    % ---------------- Callback / helper functions ----------------
    function drawScene()
        cla(hAx3); hold(hAx3,'on'); grid(hAx3,'on'); axis(hAx3,'equal');
        xlabel(hAx3,'x (mm)'); ylabel(hAx3,'y (mm)'); zlabel(hAx3,'z (mm)');
        xlim(hAx3,[0 4]); ylim(hAx3,[-4 4]); zlim(hAx3,[-6 6]);
        view(hAx3,[40 20]);

        % Nozzle (point)
        scatter3(hAx3,0,0,0,140,'k','filled'); text(hAx3,0.12,-0.25,0,'Nozzle','FontSize',10);

        % Cap plates (square lying in xy-plane at z=+/-W/2)
        x_cap_start_mm = params.cap_start_mm;
        x_cap_end_mm = x_cap_start_mm + params.L1_mm;
        y_plate_min = -params.L1_mm/2; y_plate_max = params.L1_mm/2;
        z_top_mm = params.W_mm/2; z_bot_mm = -params.W_mm/2;
        px = [x_cap_start_mm, x_cap_end_mm, x_cap_end_mm, x_cap_start_mm];
        py = [y_plate_min, y_plate_min, y_plate_max, y_plate_max];
        patch(hAx3, px, py, z_top_mm*ones(1,4), [0.8 0.9 1],'FaceAlpha',0.65,'EdgeColor','none');
        patch(hAx3, px, py, z_bot_mm*ones(1,4), [0.8 0.9 1],'FaceAlpha',0.65,'EdgeColor','none');

        % Paper at x = D_mm
        Dmm = params.D_mm;
        paper_y = [params.paper_w_mm/2, -params.paper_w_mm/2, -params.paper_w_mm/2, params.paper_w_mm/2];
        paper_z = [params.paper_h_mm/2, params.paper_h_mm/2, -params.paper_h_mm/2, -params.paper_h_mm/2];
        patch(hAx3, Dmm*ones(1,4), paper_y, paper_z, [1 1 0.9], 'FaceAlpha',0.95,'EdgeColor','k');
        scatter3(hAx3, Dmm, 0, 0, 80, 'filled', 'MarkerEdgeColor','k');
        title(hAx3,'3D View: Nozzle -> Capacitor -> Paper'); drawnow;
    end

    function computeV_cb(~,~)
        compute_and_plot_voltage();
    end

    function compute_and_plot_voltage()
        % compute Vz for I-shape and plot staircase
        mm2m = 1e-3;
        Dm = params.D_mm * mm2m;
        L1m = params.L1_mm * mm2m;
        Wm = params.W_mm * mm2m;
        cap_start_m = params.cap_start_mm * mm2m;
        cap_end_m = cap_start_m + L1m;

        r = params.droplet_d / 2;
        m = (4/3)*pi*r^3*params.rho;

        T_inside = L1m / params.vx;
        t_after = (Dm - cap_end_m) / params.vx;
        dot_spacing = (25.4e-3) / params.dpi;
        paper_h_m = params.paper_h_mm * mm2m;
        N_dots = max(1, round(paper_h_m / dot_spacing));
        z_targets = linspace(paper_h_m/2, -paper_h_m/2, N_dots);

        K = (params.q/(m*Wm)) * (0.5 * T_inside^2 + T_inside * t_after);
        V_required = z_targets / K;
        V_clip = V_required;
        V_clip(abs(V_clip) > params.Vmax) = sign(V_clip(abs(V_clip) > params.Vmax)) .* params.Vmax;

        t_fire = (0:N_dots-1) / params.f_fire;
        t_exit = t_fire + cap_end_m / params.vx;
        t_plot = [t_exit, t_exit(end) + 1/params.f_fire];
        V_plot = [V_clip, V_clip(end)];

        axes(hAxV); cla(hAxV);
        stairs(t_plot, V_plot, 'LineWidth', 1.5);
        xlabel(hAxV,'time (s)'); ylabel(hAxV,'V_z (V)'); title(hAxV,'Voltage staircase V_z(t)'); grid(hAxV,'on');

        % store sim data for run
        simData.z_targets = z_targets;
        simData.Vz = V_clip;
        simData.cap_end_m = cap_end_m;
        simData.L1m = L1m;
        simData.Wm = Wm;
        simData.N_dots = N_dots;
        simData.t_total = Dm / params.vx;
        simData.t_exit = t_exit;
        assignin('base','simData_I',simData);

        set(hStatus,'String',sprintf('Computed Vz for %d dots', N_dots));
        drawnow;
    end

    function start_cb(~,~)
        compute_and_plot_voltage();
        simData = evalin('base','simData_I');
        if isempty(simData), set(hStatus,'String','Compute V(t) first'); return; end
        if running, return; end
        running = true; stopRequested = false;
        set(hStatus,'String','Running I-print...');
        timing_metrics('start_i');                % start total I timer
        run_I_sim(simData);
        tI = timing_metrics('stop_i');            % stop total I timer; returns struct
        if ~isempty(tI) && isfield(tI,'I_time') && ~isempty(tI.I_time)
            set(hTotalTime,'String',sprintf('%.6f', tI.I_time));
        end
        running = false;
    end

    function stop_cb(~,~)
        stopRequested = true; running = false;
        set(hStatus,'String','Stop requested'); drawnow;
    end

    function reset_cb(~,~)
        stopRequested = true; running = false;
        set(hStatus,'String','Reset - ready'); set(hFirstTime,'String','N/A'); set(hTotalTime,'String','N/A');
        drawScene(); cla(hAxV);
    end

    % parameter callbacks
    function q_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.q=val; end; end
    function diam_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.droplet_d=val*1e-6; end; end
    function vx_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.vx=val; end; end
    function f_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.f_fire=val; end; end
    function vmax_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.Vmax=val; end; end
    function W_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.W_mm=val; end; end
    function capstart_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.cap_start_mm=val; end; end
    function L1_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.L1_mm=val; end; end
    function L2_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.L2_mm=val; end; end
    function D_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.D_mm=val; end; end
    function py_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.paper_w_mm=val; end; end
    function pz_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.paper_h_mm=val; end; end
    function dpi_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.dpi=val; end; end
    function speed_cb(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val); params.animation_speed=max(0.01,val); end; end

    % ---------------- Simulation core ----------------
    function run_I_sim(simData)
        mm2m = 1e-3;
        cmap = jet(simData.N_dots);
        hFlight = scatter3(hAx3, nan, nan, nan, params.markerSize_flight, 'filled', 'MarkerEdgeColor','k');
        impactHandles = gobjects(simData.N_dots,1);
        r = params.droplet_d/2; m = (4/3)*pi*r^3*params.rho;
        x_cap_start_m = params.cap_start_mm * mm2m;
        x_cap_end_m = simData.cap_end_m;

        % Loop through each droplet sequentially
        for k = 1:simData.N_dots
            if stopRequested, set(hStatus,'String','Stopped'); break; end

            % start first-drop timer at first firing
            if k == 1
                t_first_struct = timing_metrics('start_first_drop');
            end

            Vk = simData.Vz(k);
            Ez = Vk / simData.Wm;
            az_inside = (params.q * Ez) / m;

            % initial conditions
            x = 0; y = 0; z = 0;
            vx_now = params.vx; vz_now = 0;

            steps = max(30, round(100 / params.animation_speed));
            t_total = simData.t_total;
            dt_local = t_total / steps;

            for s = 1:steps
                if x >= x_cap_start_m && x <= x_cap_end_m
                    az = az_inside;
                else
                    az = 0;
                end

                vz_now = vz_now + az * dt_local;
                x = x + vx_now * dt_local;
                z = z + vz_now * dt_local;

                set(hFlight,'XData',x/mm2m,'YData',y/mm2m,'ZData',z/mm2m,'CData',cmap(k,:),'SizeData',params.markerSize_flight);
                drawnow limitrate;
            end

            % place persistent impact
            landing_z_mm = z / mm2m;
            landing_y_mm = y / mm2m;
            impactHandles(k) = scatter3(hAx3, params.D_mm, landing_y_mm, landing_z_mm, params.markerSize_impact, ...
                'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor','k');

            % stop first-drop timer when first impact occurs
            if k == 1
                t_first_result = timing_metrics('stop_first_drop');
                if ~isempty(t_first_result) && isfield(t_first_result,'first_drop') && ~isempty(t_first_result.first_drop)
                    set(hFirstTime,'String',sprintf('%.6f', t_first_result.first_drop));
                end
            end

            set(hFlight,'XData',nan,'YData',nan,'ZData',nan);
            pause(max(0, 0.002 / params.animation_speed));
        end
        set(hStatus,'String','I-print complete'); drawnow;
    end

    % ---------------- Timing utility included as nested function ----------------
    function out = timing_metrics(event)
        persistent t_first_start t_I_start
        out = struct('first_drop',[],'I_time',[]);
        switch lower(event)
            case 'start_first_drop'
                t_first_start = tic;
            case 'stop_first_drop'
                if ~isempty(t_first_start)
                    val = toc(t_first_start);
                    out.first_drop = val;
                    fprintf('First droplet time: %.6f s\n', val);
                end
            case 'start_i'
                t_I_start = tic;
            case 'stop_i'
                if ~isempty(t_I_start)
                    val = toc(t_I_start);
                    out.I_time = val;
                    fprintf('Total I draw time: %.6f s\n', val);
                end
            otherwise
                warning('Unknown timing event: %s', event);
        end
    end

end
