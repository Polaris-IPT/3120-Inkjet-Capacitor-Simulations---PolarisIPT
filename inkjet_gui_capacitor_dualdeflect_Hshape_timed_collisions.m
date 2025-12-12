function inkjet_gui_capacitor_dualdeflect_Hshape_timed_collisions
% Dual-deflection GUI (H-shape) with plate-collision detection (Option 1).
% Z-plates lie in xy-plane at z = +/- Wz/2; Y-plates lie in xz-plane at y = +/- Wy/2.
% Collisions on either plate remove the droplet and mark the collision with a red X.

    % ----------- default parameters -----------
    params.D_mm = 3.0;
    params.cap_start_mm = 1.25;
    params.L1_mm = 0.5;
    params.L2_mm = 1.25;
    params.Wz_mm = 1.0;
    params.Wy_mm = 1.0;
    params.vx = 10;
    params.droplet_d = 30e-6;
    params.rho = 1000;
    params.q = 1e-12;
    params.f_fire = 1e4;
    params.Vmax = 2000;
    params.dpi = 300;
    params.paper_w_mm = 8.5;
    params.paper_h_mm = 11.0;
    params.animation_speed = 1;
    params.markerSize_flight = 80;
    params.markerSize_impact = 60;

    running = false;
    stopRequested = false;

    % ----------- GUI -----------
    hFig = figure('Name','Inkjet Dual GUI (H) - Collisions','NumberTitle','off',...
        'Units','normalized','Position',[0.02 0.03 0.96 0.9]);

    ctrlW = 0.22; hCtrl = uipanel('Parent',hFig,'Title','Controls','FontSize',11,...
        'Units','normalized','Position',[0.01 0.02 ctrlW 0.94]);

    % 3D axes
    ax3pos = [0.26 0.12 0.5 0.82];
    hAx3 = axes('Parent',hFig,'Units','normalized','Position',ax3pos);
    hold(hAx3,'on'); grid(hAx3,'on'); axis(hAx3,'equal');
    xlabel(hAx3,'x (mm)'); ylabel(hAx3,'y (mm)'); zlabel(hAx3,'z (mm)');
    view(hAx3,[40 20]); xlim(hAx3,[0 4]); ylim(hAx3,[-4 4]); zlim(hAx3,[-6 6]);

    % Voltage axes stacked
    axVposTop = [0.78 0.52 0.20 0.42];
    axVposBot = [0.78 0.12 0.20 0.36];
    hAxVz = axes('Parent',hFig,'Units','normalized','Position',axVposTop);
    title(hAxVz,'Vz(t)'); xlabel(hAxVz,'time (s)'); ylabel(hAxVz,'V (V)'); grid(hAxVz,'on');
    hAxVy = axes('Parent',hFig,'Units','normalized','Position',axVposBot);
    title(hAxVy,'Vy(t)'); xlabel(hAxVy,'time (s)'); ylabel(hAxVy,'V (V)'); grid(hAxVy,'on');

    % Controls layout
    ctrlY = 0.92; step = 0.055;
    addLabel = @(txt,y) uicontrol('Parent',hCtrl,'Style','text','String',txt,...
        'Units','normalized','Position',[0.02 y 0.46 0.04],'HorizontalAlignment','left');
    addEdit = @(val,y,cb) uicontrol('Parent',hCtrl,'Style','edit','String',num2str(val),...
        'Units','normalized','Position',[0.5 y 0.46 0.05],'Callback',cb);

    % NOTE: callbacks use unique names cb_* to avoid duplicate nested function names
    addLabel('Droplet charge q (C):', ctrlY); hQ = addEdit(params.q, ctrlY, @cb_q); ctrlY = ctrlY - step;
    addLabel('Droplet diameter (um):', ctrlY); hDiam = addEdit(params.droplet_d*1e6, ctrlY, @cb_diam); ctrlY = ctrlY - step;
    addLabel('Droplet speed vx (m/s):', ctrlY); hVx = addEdit(params.vx, ctrlY, @cb_vx); ctrlY = ctrlY - step;
    addLabel('Firing freq f (Hz):', ctrlY); hF = addEdit(params.f_fire, ctrlY, @cb_f); ctrlY = ctrlY - step;
    addLabel('Max voltage Vmax (V):', ctrlY); hVmax = addEdit(params.Vmax, ctrlY, @cb_vmax); ctrlY = ctrlY - step;
    addLabel('Z-plate separation Wz (mm):', ctrlY); hWz = addEdit(params.Wz_mm, ctrlY, @cb_Wz); ctrlY = ctrlY - step;
    addLabel('Y-plate separation Wy (mm):', ctrlY); hWy = addEdit(params.Wy_mm, ctrlY, @cb_Wy); ctrlY = ctrlY - step;
    addLabel('Cap start x (mm):', ctrlY); hCapStart = addEdit(params.cap_start_mm, ctrlY, @cb_capstart); ctrlY = ctrlY - step;
    addLabel('L1 (mm):', ctrlY); hL1 = addEdit(params.L1_mm, ctrlY, @cb_L1); ctrlY = ctrlY - step;
    addLabel('L2 (mm):', ctrlY); hL2 = addEdit(params.L2_mm, ctrlY, @cb_L2); ctrlY = ctrlY - step;
    addLabel('Paper x (mm):', ctrlY); hD = addEdit(params.D_mm, ctrlY, @cb_D); ctrlY = ctrlY - step;
    addLabel('Paper width (y) mm:', ctrlY); hPy = addEdit(params.paper_w_mm, ctrlY, @cb_py); ctrlY = ctrlY - step;
    addLabel('Paper height (z) mm:', ctrlY); hPz = addEdit(params.paper_h_mm, ctrlY, @cb_pz); ctrlY = ctrlY - step;
    addLabel('DPI:', ctrlY); hDPI = addEdit(params.dpi, ctrlY, @cb_dpi); ctrlY = ctrlY - step;
    addLabel('Animation speed (mult):', ctrlY); hSpeed = addEdit(params.animation_speed, ctrlY, @cb_speed); ctrlY = ctrlY - step;

    % Buttons & status/timing
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Start H-print','Units','normalized',...
        'Position',[0.02 0.05 0.22 0.06],'BackgroundColor',[0.6 1 0.6],'Callback',@start_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Stop','Units','normalized',...
        'Position',[0.27 0.05 0.22 0.06],'BackgroundColor',[1 0.6 0.6],'Callback',@stop_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Reset','Units','normalized',...
        'Position',[0.52 0.05 0.22 0.06],'Callback',@reset_cb);
    uicontrol('Parent',hCtrl,'Style','pushbutton','String','Compute V(t)','Units','normalized',...
        'Position',[0.02 0.01 0.96 0.035],'Callback',@computeV_cb);

    hStatus = uicontrol('Parent',hCtrl,'Style','text','String','Ready','Units','normalized',...
        'Position',[0.02 0.12 0.96 0.045],'BackgroundColor',get(hFig,'Color'),'HorizontalAlignment','left');

    addLabel('First droplet time (s):', 0.18);
    hFirstTime = uicontrol('Parent',hCtrl,'Style','text','String','N/A','Units','normalized',...
        'Position',[0.5 0.18 0.46 0.04],'HorizontalAlignment','left');
    addLabel('Total H draw time (s):', 0.13);
    hTotalTime = uicontrol('Parent',hCtrl,'Style','text','String','N/A','Units','normalized',...
        'Position',[0.5 0.13 0.46 0.04],'HorizontalAlignment','left');

    drawScene();

    % ---------------- scene draw ----------------
    function drawScene()
        cla(hAx3); hold(hAx3,'on'); grid(hAx3,'on'); axis(hAx3,'equal');
        xlabel(hAx3,'x (mm)'); ylabel(hAx3,'y (mm)'); zlabel(hAx3,'z (mm)');
        xlim(hAx3,[0 4]); ylim(hAx3,[-4 4]); zlim(hAx3,[-6 6]);
        view(hAx3,[40 20]);

        scatter3(hAx3,0,0,0,140,'k','filled'); text(hAx3,0.12,-0.25,0,'Nozzle','FontSize',10);

        x_cap_start_mm = params.cap_start_mm; x_cap_end_mm = x_cap_start_mm + params.L1_mm;
        y_plate_min = -params.L1_mm/2; y_plate_max = params.L1_mm/2;
        z_top_mm = params.Wz_mm/2; z_bot_mm = -params.Wz_mm/2;
        px = [x_cap_start_mm, x_cap_end_mm, x_cap_end_mm, x_cap_start_mm];
        py = [y_plate_min, y_plate_min, y_plate_max, y_plate_max];

        % Z-plates
        patch(hAx3, px, py, z_top_mm*ones(1,4), [0.8 0.9 1],'FaceAlpha',0.65,'EdgeColor','none');
        patch(hAx3, px, py, z_bot_mm*ones(1,4), [0.8 0.9 1],'FaceAlpha',0.65,'EdgeColor','none');

        % Y-plates (xz-plane)
        z_pl_min = -params.L1_mm/2; z_pl_max = params.L1_mm/2;
        pz = [z_pl_min, z_pl_min, z_pl_max, z_pl_max];
        patch(hAx3, px, params.Wy_mm/2*ones(1,4), pz, [0.9 0.8 1], 'FaceAlpha',0.65,'EdgeColor','none');
        patch(hAx3, px, -params.Wy_mm/2*ones(1,4), pz, [0.9 0.8 1], 'FaceAlpha',0.65,'EdgeColor','none');

        % paper at x = D
        Dmm = params.D_mm;
        paper_y = [params.paper_w_mm/2, -params.paper_w_mm/2, -params.paper_w_mm/2, params.paper_w_mm/2];
        paper_z = [params.paper_h_mm/2, params.paper_h_mm/2, -params.paper_h_mm/2, -params.paper_h_mm/2];
        patch(hAx3, Dmm*ones(1,4), paper_y, paper_z, [1 1 0.9], 'FaceAlpha',0.95,'EdgeColor','k');
        scatter3(hAx3, Dmm, 0, 0, 80, 'filled', 'MarkerEdgeColor','k');
        title(hAx3,'3D View: Nozzle -> Dual plates -> Paper'); drawnow;
    end

    % ---------------- compute voltages & staircases ----------------
    function computeV_cb(~,~)
        mm2m = 1e-3;
        Dm = params.D_mm * mm2m;
        L1m = params.L1_mm * mm2m;
        Wzm = params.Wz_mm * mm2m;
        Wym = params.Wy_mm * mm2m;
        cap_start_m = params.cap_start_mm * mm2m;
        cap_end_m = cap_start_m + L1m;

        r = params.droplet_d / 2;
        m = (4/3) * pi * r^3 * params.rho;

        T_inside = L1m / params.vx;
        t_after = (Dm - cap_end_m) / params.vx;

        dot_spacing_m = (25.4e-3) / params.dpi;
        paper_h_m = params.paper_h_mm * mm2m;

        % build H path
        y_left = -4e-3;
        z_vert = linspace(paper_h_m/2, -paper_h_m/2, max(2, round(paper_h_m/dot_spacing_m)));
        left_pts = [repmat(y_left, size(z_vert)) ; z_vert]';
        z_mid = 0;
        y_horiz = linspace(-4e-3, 4e-3, max(2, round((8e-3)/dot_spacing_m)));
        mid_pts = [y_horiz ; repmat(z_mid, size(y_horiz))]';
        y_right = 4e-3;
        right_pts = [repmat(y_right, size(z_vert)); z_vert]';
        targets_yz = [left_pts; mid_pts; right_pts];
        N_dots = size(targets_yz,1);

        Kz = (params.q / (m * Wzm)) * (0.5 * T_inside^2 + T_inside * t_after);
        Ky = (params.q / (m * Wym)) * (0.5 * T_inside^2 + T_inside * t_after);

        z_targets = targets_yz(:,2); y_targets = targets_yz(:,1);
        Vz_req = z_targets / Kz; Vy_req = y_targets / Ky;
        Vz_clip = Vz_req; Vy_clip = Vy_req;
        Vz_clip(abs(Vz_clip) > params.Vmax) = sign(Vz_clip(abs(Vz_clip)>params.Vmax)) .* params.Vmax;
        Vy_clip(abs(Vy_clip) > params.Vmax) = sign(Vy_clip(abs(Vy_clip)>params.Vmax)) .* params.Vmax;

        t_fire = (0:N_dots-1) / params.f_fire;
        t_exit = t_fire + cap_end_m / params.vx;
        t_plot = [t_exit, t_exit(end) + 1/params.f_fire];

        axes(hAxVz); cla(hAxVz);
        stairs(t_plot, [Vz_clip; Vz_clip(end)], 'LineWidth', 1.4);
        xlabel(hAxVz,'time (s)'); ylabel(hAxVz,'Vz (V)'); title(hAxVz,'Vz(t)'); grid(hAxVz,'on');

        axes(hAxVy); cla(hAxVy);
        stairs(t_plot, [Vy_clip; Vy_clip(end)], 'LineWidth', 1.4);
        xlabel(hAxVy,'time (s)'); ylabel(hAxVy,'Vy (V)'); title(hAxVy,'Vy(t)'); grid(hAxVy,'on');

        simData.targets_yz = targets_yz;
        simData.Vz = Vz_clip; simData.Vy = Vy_clip;
        simData.cap_start_m = cap_start_m; simData.cap_end_m = cap_end_m;
        simData.L1m = L1m; simData.Wzm = Wzm; simData.Wym = Wym;
        simData.N_dots = N_dots; simData.t_total = Dm / params.vx;
        assignin('base','simData_H',simData);

        set(hStatus,'String',sprintf('Computed voltages for %d dots', N_dots)); drawnow;
    end

    % ---------------- start / stop / reset ----------------
    function start_cb(~,~)
        computeV_cb();
        simData = evalin('base','simData_H');
        if isempty(simData), set(hStatus,'String','Compute V(t) first'); return; end
        if running, return; end
        running = true; stopRequested = false;
        set(hStatus,'String','Running H-print...');
        timing_metrics('start_h');
        run_H_sim(simData);
        tH = timing_metrics('stop_h');
        if ~isempty(tH) && isfield(tH,'H_time') && ~isempty(tH.H_time)
            set(hTotalTime,'String',sprintf('%.6f', tH.H_time));
        end
        running = false;
    end

    function stop_cb(~,~)
        stopRequested = true; running = false; set(hStatus,'String','Stop requested'); drawnow;
    end

    function reset_cb(~,~)
        stopRequested = true; running = false;
        set(hStatus,'String','Reset - ready'); set(hFirstTime,'String','N/A'); set(hTotalTime,'String','N/A');
        drawScene(); cla(hAxVz); cla(hAxVy); drawnow;
    end

    % ---------------- param callbacks (renamed cb_*) ----------------
    function cb_q(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.q=val; end; end
    function cb_diam(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.droplet_d=val*1e-6; end; end
    function cb_vx(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.vx=val; end; end
    function cb_f(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.f_fire=val; end; end
    function cb_vmax(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.Vmax=val; end; end
    function cb_Wz(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.Wz_mm=val; end; end
    function cb_Wy(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.Wy_mm=val; end; end
    function cb_capstart(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.cap_start_mm=val; end; end
    function cb_L1(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.L1_mm=val; end; end
    function cb_L2(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.L2_mm=val; end; end
    function cb_D(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.D_mm=val; end; end
    function cb_py(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.paper_w_mm=val; end; end
    function cb_pz(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.paper_h_mm=val; end; end
    function cb_dpi(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.dpi=val; end; end
    function cb_speed(hObj,~); val=str2double(get(hObj,'String')); if ~isnan(val), params.animation_speed=max(0.01,val); end; end

    % ---------------- core H simulation (with collisions) ----------------
    function run_H_sim(simData)
        mm2m = 1e-3;
        cmap = jet(simData.N_dots);
        hFlight = scatter3(hAx3, nan, nan, nan, params.markerSize_flight, 'filled', 'MarkerEdgeColor','k');
        impactHandles = gobjects(simData.N_dots,1);
        collisionHandles = gobjects(simData.N_dots,1);

        r_drop = params.droplet_d / 2;
        m = (4/3) * pi * r_drop^3 * params.rho;

        x_cap_start_m = simData.cap_start_m;
        x_cap_end_m = simData.cap_end_m;

        % plate footprints
        y_plate_min_m = - (params.L1_mm/2) * mm2m;
        y_plate_max_m =   (params.L1_mm/2) * mm2m;
        z_plate_min_m = - (params.L1_mm/2) * mm2m;
        z_plate_max_m =   (params.L1_mm/2) * mm2m;

        z_top_m = (params.Wz_mm/2) * mm2m;
        z_bot_m = -(params.Wz_mm/2) * mm2m;
        y_top_m = (params.Wy_mm/2) * mm2m;
        y_bot_m = -(params.Wy_mm/2) * mm2m;

        first_drop_recorded = false;

        for k = 1:simData.N_dots
            if stopRequested, set(hStatus,'String','Stopped'); break; end

            % start-first-drop timer at first firing
            if k == 1
                timing_metrics('start_first_drop');
            end

            Vz = simData.Vz(k);
            Vy = simData.Vy(k);

            az_inside = (params.q * (Vz / simData.Wzm)) / m;
            ay_inside = (params.q * (Vy / simData.Wym)) / m;

            x = 0; y = 0; z = 0;
            vx_now = params.vx; vy_now = 0; vz_now = 0;

            steps = max(40, round(120 / params.animation_speed));
            t_total = simData.t_total;
            dt_local = t_total / steps;
            collided = false;

            for s = 1:steps
                if x >= x_cap_start_m && x <= x_cap_end_m
                    az = az_inside; ay = ay_inside;
                else
                    az = 0; ay = 0;
                end

                vz_now = vz_now + az * dt_local;
                vy_now = vy_now + ay * dt_local;
                x = x + vx_now * dt_local;
                z = z + vz_now * dt_local;
                y = y + vy_now * dt_local;

                % collision check: Z-plates (top/bottom)
                if (x >= x_cap_start_m && x <= x_cap_end_m) && (y >= y_plate_min_m && y <= y_plate_max_m)
                    if z + r_drop >= z_top_m
                        collided = true;
                        coll_x_mm = x / mm2m; coll_y_mm = y / mm2m; coll_z_mm = z_top_m / mm2m;
                        collisionHandles(k) = scatter3(hAx3, coll_x_mm, coll_y_mm, coll_z_mm, 80, 'x', 'MarkerEdgeColor','r','LineWidth',2);
                        break;
                    end
                    if z - r_drop <= z_bot_m
                        collided = true;
                        coll_x_mm = x / mm2m; coll_y_mm = y / mm2m; coll_z_mm = z_bot_m / mm2m;
                        collisionHandles(k) = scatter3(hAx3, coll_x_mm, coll_y_mm, coll_z_mm, 80, 'x', 'MarkerEdgeColor','r','LineWidth',2);
                        break;
                    end
                end

                % collision check: Y-plates (left/right), plates lie at y = +/- Wy/2
                if (x >= x_cap_start_m && x <= x_cap_end_m) && (z >= z_plate_min_m && z <= z_plate_max_m)
                    if y + r_drop >= y_top_m
                        collided = true;
                        coll_x_mm = x / mm2m; coll_y_mm = y_top_m / mm2m; coll_z_mm = z / mm2m;
                        collisionHandles(k) = scatter3(hAx3, coll_x_mm, coll_y_mm, coll_z_mm, 80, 'x', 'MarkerEdgeColor','r','LineWidth',2);
                        break;
                    end
                    if y - r_drop <= y_bot_m
                        collided = true;
                        coll_x_mm = x / mm2m; coll_y_mm = y_bot_m / mm2m; coll_z_mm = z / mm2m;
                        collisionHandles(k) = scatter3(hAx3, coll_x_mm, coll_y_mm, coll_z_mm, 80, 'x', 'MarkerEdgeColor','r','LineWidth',2);
                        break;
                    end
                end

                set(hFlight,'XData',x/mm2m,'YData',y/mm2m,'ZData',z/mm2m,'CData',cmap(k,:),'SizeData',params.markerSize_flight);
                drawnow limitrate;
            end

            if collided
                set(hFlight,'XData',nan,'YData',nan,'ZData',nan);
                pause(max(0, 0.002 / params.animation_speed));
                continue;
            end

            % successful impact: place persistent marker
            landing_y_mm = y / mm2m; landing_z_mm = z / mm2m;
            impactHandles(k) = scatter3(hAx3, params.D_mm, landing_y_mm, landing_z_mm, params.markerSize_impact, ...
                'MarkerFaceColor', cmap(k,:), 'MarkerEdgeColor','k');

            if ~first_drop_recorded
                t_first_result = timing_metrics('stop_first_drop');
                if ~isempty(t_first_result) && isfield(t_first_result,'first_drop') && ~isempty(t_first_result.first_drop)
                    set(hFirstTime,'String',sprintf('%.6f', t_first_result.first_drop));
                end
                first_drop_recorded = true;
            end

            set(hFlight,'XData',nan,'YData',nan,'ZData',nan);
            pause(max(0, 0.002 / params.animation_speed));
        end

        set(hStatus,'String','H-print complete'); drawnow;
    end

    % ---------------- timing utility nested function ----------------
    function out = timing_metrics(event)
        persistent t_first_start t_H_start
        out = struct('first_drop',[],'H_time',[]);
        switch lower(event)
            case 'start_first_drop'
                t_first_start = tic;
            case 'stop_first_drop'
                if ~isempty(t_first_start)
                    val = toc(t_first_start);
                    out.first_drop = val;
                    fprintf('First droplet time: %.6f s\n', val);
                end
            case 'start_h'
                t_H_start = tic;
            case 'stop_h'
                if ~isempty(t_H_start)
                    val = toc(t_H_start);
                    out.H_time = val;
                    fprintf('Total H draw time: %.6f s\n', val);
                end
            otherwise
                warning('Unknown timing event: %s', event);
        end
    end

end
