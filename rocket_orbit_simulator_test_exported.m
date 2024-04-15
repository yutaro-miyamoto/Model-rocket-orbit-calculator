classdef rocket_orbit_simulator_test_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        start_button                    matlab.ui.control.StateButton
        Panel_launch_deg                matlab.ui.container.Panel
        phi                             matlab.ui.control.NumericEditField
        LaunchdirectionphidegLabel      matlab.ui.control.Label
        Image                           matlab.ui.control.Image
        theta                           matlab.ui.control.NumericEditField
        ZdegLabel                       matlab.ui.control.Label
        InitializeallgraphsButton       matlab.ui.control.Button
        simulation_details              matlab.ui.container.Panel
        num_simulations                 matlab.ui.control.NumericEditField
        Label                           matlab.ui.control.Label
        simulation_time                 matlab.ui.control.NumericEditField
        SimulationtimesLabel            matlab.ui.control.Label
        fundamental_details             matlab.ui.container.Panel
        wind_effectness                 matlab.ui.control.NumericEditField
        Label_10                        matlab.ui.control.Label
        wind_degree                     matlab.ui.control.NumericEditField
        WinddirectiondegLabel           matlab.ui.control.Label
        wind_speed                      matlab.ui.control.NumericEditField
        Airvelocityms1Label             matlab.ui.control.Label
        gravitational_acceleration      matlab.ui.control.NumericEditField
        GravitationalaccelerationLabel  matlab.ui.control.Label
        rho                             matlab.ui.control.NumericEditField
        Label_3                         matlab.ui.control.Label
        Panel_parachute                 matlab.ui.container.Panel
        C_xo                            matlab.ui.control.NumericEditField
        Label_8                         matlab.ui.control.Label
        drag_coefficient_parachute      matlab.ui.control.NumericEditField
        Label_7                         matlab.ui.control.Label
        parachute_S                     matlab.ui.control.NumericEditField
        Label_6                         matlab.ui.control.Label
        parachute_time                  matlab.ui.control.NumericEditField
        Label_5                         matlab.ui.control.Label
        parachute_on_off                matlab.ui.control.DropDown
        figures                         matlab.ui.container.TabGroup
        ThefinaldestinationTab          matlab.ui.container.Tab
        final_destination               matlab.ui.control.UIAxes
        OrbitTab                        matlab.ui.container.Tab
        orbit                           matlab.ui.control.UIAxes
        vtgraphTab                      matlab.ui.container.Tab
        v_t_graph                       matlab.ui.control.UIAxes
        htgraphTab                      matlab.ui.container.Tab
        h_t_graph                       matlab.ui.control.UIAxes
        vhtTab                          matlab.ui.container.Tab
        v_h_t_graph                     matlab.ui.control.UIAxes
        rocket                          matlab.ui.container.Panel
        WeightunitsettingLabel          matlab.ui.control.Label
        cross_sectional_side_area       matlab.ui.control.NumericEditField
        Label_12                        matlab.ui.control.Label
        fuel_weight_end                 matlab.ui.control.NumericEditField
        kgLabel                         matlab.ui.control.Label
        total_impulse                   matlab.ui.control.NumericEditField
        NsecLabel                       matlab.ui.control.Label
        fuel_weight_begin               matlab.ui.control.NumericEditField
        Label_11                        matlab.ui.control.Label
        thrust_time                     matlab.ui.control.NumericEditField
        InjectiontimesLabel             matlab.ui.control.Label
        drag_coefficient                matlab.ui.control.NumericEditField
        Label_2                         matlab.ui.control.Label
        constant_weight                 matlab.ui.control.NumericEditField
        gLabel                          matlab.ui.control.Label
        cross_sectional_area            matlab.ui.control.NumericEditField
        Label_4                         matlab.ui.control.Label
        checking_kg_or_g                matlab.ui.control.DropDown
    end


    % Public properties that correspond to the Simulink model
    properties (Access = public, Transient)
        Simulation simulink.Simulation
    end

 

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % 外部関数があるディレクトリを追加
            addpath('/class_dir');
        end

        % Value changed function: num_simulations
        function num_simulationsValueChanged(app, event)
            %シミュレーション回数
            value = app.num_simulations.Value;            
        end

        % Value changed function: gravitational_acceleration
        function gravitational_accelerationValueChanged(app, event)
            %重力加速度
            value = app.gravitational_acceleration.Value;            
        end

        % Value changed function: constant_weight
        function constant_weightValueChanged(app, event)
            % 固定質量/kg
            value = app.constant_weight.Value;
        end

        % Value changed function: fuel_weight_end
        function fuel_weight_endValueChanged2(app, event)
            %燃焼開始時燃料質量/kg
            value = app.fuel_weight_end.Value;            
        end

        % Value changed function: fuel_weight_begin
        function fuel_weight_beginValueChanged(app, event)
            % 燃焼終了時質量
            value = app.fuel_weight_begin.Value;           
        end

        % Value changed function: cross_sectional_area
        function cross_sectional_areaValueChanged(app, event)
            % 機体断面積
            value = app.cross_sectional_area.Value;            
        end

        % Value changed function: drag_coefficient
        function drag_coefficientValueChanged(app, event)
            % 機体の抵抗係数
            value = app.drag_coefficient.Value;            
        end

        % Value changed function: rho
        function rhoValueChanged(app, event)
            % 空気密度
            value = app.rho.Value;
        end

        % Value changed function: total_impulse
        function total_impulseValueChanged(app, event)
            %トータルインパルス
            value = app.total_impulse.Value;            
        end

        % Value changed function: thrust_time
        function thrust_timeValueChanged(app, event)
            %噴射時間
            value = app.thrust_time.Value;            
        end

        % Value changed function: theta
        function thetaValueChanged(app, event)
            % Z軸とベクトルのなす角(0<=theta<=180)
            value = app.theta.Value;            
        end

        % Value changed function: phi
        function phiValueChanged(app, event)
            % 打ち上げ角度
            % X軸とxy座標上でのベクトルの正射影となす角(0<=phi<=360)
            value = app.phi.Value;            
        end

        % Value changed function: simulation_time
        function simulation_timeValueChanged(app, event)
            % シミュレーション時間/s
            value = app.simulation_time.Value;            
        end

        % Value changed function: wind_speed
        function wind_speedValueChanged(app, event)
            % 風速度
            value = app.wind_speed.Value;            
        end

        % Value changed function: wind_degree
        function wind_degreeValueChanged(app, event)
            %風向き入力
            value = app.wind_degree.Value;        
        end

        % Value changed function: wind_effectness
        function wind_effectnessValueChanged(app, event)
            % 風の影響度
            value = app.wind_effectness.Value;            
        end

        % Value changed function: parachute_time
        function parachute_timeValueChanged(app, event)
            %パラシュート展開時間
            value = app.parachute_time.Value;            
        end

        % Value changed function: parachute_S
        function parachute_SValueChanged(app, event)
            % パラシュートの投影面積
            value = app.parachute_S.Value;            
        end

        % Value changed function: drag_coefficient_parachute
        function drag_coefficient_parachuteValueChanged(app, event)
            % パラシュートのCd
            value = app.drag_coefficient_parachute.Value;            
        end

        % Value changed function: C_xo
        function C_xoValueChanged(app, event)
            %傘体荷重係数の導入
            value = app.C_xo.Value;            
        end

        % Value changed function: checking_kg_or_g
        function checking_kg_or_gValueChanged(app, event)
            value = app.checking_kg_or_g.Value;
        end

        % Value changed function: parachute_on_off
        function parachute_on_offValueChanged(app, event)
            value = app.parachute_on_off.Value;
        end

        % Value changed function: cross_sectional_side_area
        function cross_sectional_side_areaValueChanged(app, event)
            value = app.cross_sectional_side_area.Value;
        end

        % Button pushed function: InitializeallgraphsButton
        function InitializeallgraphsButtonPushed(app, event)
            cla(app.orbit); % Axes コンポーネントをクリア
            cla(app.final_destination);
            cla(app.v_t_graph);
            cla(app.h_t_graph);
            cla(app.v_h_t_graph);
            yyaxis(app.v_h_t_graph,"left");
            cla(app.v_h_t_graph); % 左の軸を初期化
        end

        % Callback function: orbit, start_button
        function ButtonValueChanged(app, event)
            try
                % シミュレーションの時間/s
                tmax = app.simulation_time.Value;

                %kg or g の決定
                if strcmpi(app.checking_kg_or_g.Value, 'kg')
                    times_const = 1;
                elseif strcmpi(app.checking_kg_or_g.Value, 'g')
                    times_const = 0.001;
                end
                %パラシュートのon/offを切り替え
                if strcmpi(app.parachute_on_off.Value,'Parachute off')
                    app.parachute_time.Value = tmax;
                end

                % シミュレーションの開始
                num_simulations = app.num_simulations.Value;

                % 最終到着地点を保存する配列
                landing_points = zeros(2, num_simulations);

                for a = 1:num_simulations
                    % 初期位置
                    r_0 = [0;0;0];
                    r = r_0;

                    % 時間ステップ/s
                    dt = 0.001;

                    % various_setting関数の確認
                    setting = various_setting(app.fuel_weight_begin.Value,app.fuel_weight_end.Value,times_const,app.constant_weight.Value,app.thrust_time.Value,app.total_impulse.Value,dt,app.gravitational_acceleration.Value);
                    m = setting.weight;

                    % 初期速度を打ち上げ角度に基づいて計算
                    v_0 = [0;0;0];
                    v = v_0;

                    heading_vector = [sind(app.theta.Value)*cosd(app.phi.Value);sind(app.theta.Value)*sind(app.phi.Value);cosd(app.theta.Value)];
                    % setting.F_rv = F_r0(3) * heading_vector;

                    % 軌道データを保存する配列
                    r_data = zeros(3, tmax / dt);

                    % 初期位置を保存
                    r_data(:,1) = r_0;

                    % 時間データを保存する配列
                    t_data = zeros(1, tmax / dt);
                    t_data(1) = 0;

                    % 速度データを保存する配列
                    v_data = zeros(1, tmax / dt);
                    v_data(1) = 0;
                    
                    %空気抵抗
                    Air_resistance = [0;0;0];

                    % ループを実行
                    % 風速度
                    % wind_speed = [app.wind_speed.Value*cosd(app.wind_degree.Value) + app.wind_effectness.Value * rand;...
                    % app.wind_speed.Value*sind(app.wind_degree.Value)+ app.wind_effectness.Value *rand;0];
                   
                    % wind_speed = [0;0;0];%テスト用

                    for t = 1:tmax/dt   %注意！！！この中の計算はすべてdt秒あたり、で考えないと値がおかしくなる。
                        % % 風速度
                        % wind_speed = [app.wind_speed.Value*cosd(app.wind_degree.Value);...
                        %     app.wind_speed.Value*sind(app.wind_degree.Value);0];
                        

                        if t <= app.thrust_time.Value / dt %%ロケットついてるとき
                            
                            wind_speed = [app.wind_speed.Value*cosd(app.wind_degree.Value) * app.wind_effectness.Value * rand;...
                                app.wind_speed.Value*sind(app.wind_degree.Value)* app.wind_effectness.Value *rand;0];
                            %力学計算
                            if r(3)<1
                                wind_speed = [0;0;0];
                                elements = calculation(m,app.gravitational_acceleration.Value,Air_resistance,v_0,r_0,setting.weight_fuel_using,wind_speed,dt,setting.Isp,heading_vector,app.drag_coefficient.Value,app.rho.Value,app.cross_sectional_area.Value,app.cross_sectional_side_area.Value);
                                v = elements.velocity;
                                r = elements.location;

                                % 質量の更新(燃料の消費)
                                m1 = m;
                                m = m1 - setting.weight_fuel_using * dt;%dtあたりの燃料消費量を引く
                            else

                                elements = calculation(m,app.gravitational_acceleration.Value,Air_resistance,v_0,r_0,setting.weight_fuel_using,wind_speed,dt,setting.Isp,heading_vector,app.drag_coefficient.Value,app.rho.Value,app.cross_sectional_area.Value,app.cross_sectional_side_area.Value);
                                v = elements.velocity;
                                r = elements.location;

                                % 質量の更新(燃料の消費)
                                m1 = m;
                                m = m1 - setting.weight_fuel_using * dt;%dtあたりの燃料消費量を引く
                            end
                        
                        elseif  t > app.thrust_time.Value / dt && t <= (app.thrust_time.Value + app.parachute_time.Value) / dt
                            %%ロケット推進消えたE-
                            %力学計算
                            elements = calculation(m,app.gravitational_acceleration.Value,Air_resistance,v_0,r_0,0,wind_speed,dt,0,heading_vector,app.drag_coefficient.Value,app.rho.Value,app.cross_sectional_area.Value,app.cross_sectional_side_area.Value);
                            v = elements.velocity;
                            r = elements.location;

                        else%パラシュート展開。抗力係数と機体の投影面積が変わる。
                            %空気抵抗の更新
                            %https://nociws.github.io/parachute/
                            %を参照
                            Air_resistance = app.drag_coefficient_parachute.Value * 0.5 * app.rho.Value * sqrt(sum(v.^2))...
                                * v * app.parachute_S.Value * app.C_xo.Value;

                            %力学計算
                            elements = calculation(m,app.gravitational_acceleration.Value,Air_resistance,v_0,r_0,0,wind_speed,dt,0,heading_vector,app.drag_coefficient.Value,app.rho.Value,app.cross_sectional_area.Value,app.cross_sectional_side_area.Value);
                            v = elements.velocity;
                            r = elements.location;
                        end

                        % rのz座標がゼロになった場合、ループから抜ける
                        if r(3)<0
                            % z座標が無視されることを保証するために、rの3番目の要素を0にする
                            r(3) = 0;
                            % 結果を保存
                            landing_points(:,a) = r(1:2);% x座標とy座標のみを保存する
                            break
                        end


                        % 次の時間ステップの準備

                        heading_vector = elements.heading_vector;
                        Air_resistance = elements.air_resistance;

                        % 結果を保存
                        r_data(:,t+1) = r;

                        v_data(t+1) = elements.mag;

                        t_data(t+1) = t*dt;

                        % F_rv = F_r0(3) * heading_vector;%動かなくなるよ～

                        %Dv = app.total_impulse.Value / m_fuel * heading_vector;

                        r_0= r;

                        v_0 = v;

                    end


                    plot3(app.orbit, r_data(1,1:t), r_data(2,1:t), r_data(3,1:t));
                    hold(app.orbit, "on"); % Axes コンポーネントに hold on を設定

                    plot(app.v_t_graph, t_data(1,1:t), v_data(1,1:t),"-");
                    hold(app.v_t_graph, "on");

                    plot(app.h_t_graph, t_data(1,1:t), r_data(3,1:t),"-");
                    hold(app.h_t_graph, "on");

                    yyaxis(app.v_h_t_graph, 'left');
                    plot(app.v_h_t_graph, t_data(1,1:t), v_data(1,1:t), '-');
                    hold(app.v_h_t_graph, "on");
                    ylim(app.v_h_t_graph, [0, max(v_data(1,1:t))]);
                    ylabel(app.v_h_t_graph, 'Velocity/ms^{-1}');

                    yyaxis(app.v_h_t_graph, 'right');
                    plot(app.v_h_t_graph, t_data(1,1:t), r_data(3,1:t), ':');
                    hold(app.v_h_t_graph, "on");
                    ylim(app.v_h_t_graph, [0, max(r_data(3,1:t))]);
                    ylabel(app.v_h_t_graph, 'Height/m');

                end
                % 着陸地点をプロットするコード
                scatter(app.final_destination, landing_points(1,:), landing_points(2,:), 'red','*');
                hold(app.final_destination,"on");

            catch ME
                % エラーが発生した場合の処理
                errordlg(['Error starting simulation: ' ME.message], 'Simulation Error');
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Get the file path for locating images
            pathToMLAPP = fileparts(mfilename('fullpath'));

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 963 629];
            app.UIFigure.Name = 'MATLAB App';

            % Create rocket
            app.rocket = uipanel(app.UIFigure);
            app.rocket.Title = 'Configuration';
            app.rocket.Position = [10 445 473 174];

            % Create checking_kg_or_g
            app.checking_kg_or_g = uidropdown(app.rocket);
            app.checking_kg_or_g.Items = {'kg', 'g'};
            app.checking_kg_or_g.ValueChangedFcn = createCallbackFcn(app, @checking_kg_or_gValueChanged, true);
            app.checking_kg_or_g.FontSize = 8;
            app.checking_kg_or_g.Position = [395 10 47 22];
            app.checking_kg_or_g.Value = 'g';

            % Create Label_4
            app.Label_4 = uilabel(app.rocket);
            app.Label_4.HorizontalAlignment = 'right';
            app.Label_4.FontSize = 8;
            app.Label_4.Interpreter = 'tex';
            app.Label_4.Position = [-2 116 190 22];
            app.Label_4.Text = 'Area of longitudinal section of aircraft/m^{2}';

            % Create cross_sectional_area
            app.cross_sectional_area = uieditfield(app.rocket, 'numeric');
            app.cross_sectional_area.ValueChangedFcn = createCallbackFcn(app, @cross_sectional_areaValueChanged, true);
            app.cross_sectional_area.FontSize = 8;
            app.cross_sectional_area.Position = [196 114 52 28];
            app.cross_sectional_area.Value = 0.000625;

            % Create gLabel
            app.gLabel = uilabel(app.rocket);
            app.gLabel.HorizontalAlignment = 'right';
            app.gLabel.FontSize = 8;
            app.gLabel.Interpreter = 'tex';
            app.gLabel.Position = [301 89 69 22];
            app.gLabel.Text = 'Constant Mass';

            % Create constant_weight
            app.constant_weight = uieditfield(app.rocket, 'numeric');
            app.constant_weight.ValueChangedFcn = createCallbackFcn(app, @constant_weightValueChanged, true);
            app.constant_weight.FontSize = 8;
            app.constant_weight.Position = [399 87 38 25];
            app.constant_weight.Value = 26.3;

            % Create Label_2
            app.Label_2 = uilabel(app.rocket);
            app.Label_2.HorizontalAlignment = 'right';
            app.Label_2.FontSize = 8;
            app.Label_2.Interpreter = 'tex';
            app.Label_2.Position = [18 71 165 34];
            app.Label_2.Text = 'Coefficient of fuselage resistance(Cd)';

            % Create drag_coefficient
            app.drag_coefficient = uieditfield(app.rocket, 'numeric');
            app.drag_coefficient.ValueChangedFcn = createCallbackFcn(app, @drag_coefficientValueChanged, true);
            app.drag_coefficient.FontSize = 8;
            app.drag_coefficient.Position = [209 78 39 27];
            app.drag_coefficient.Value = 0.3;

            % Create InjectiontimesLabel
            app.InjectiontimesLabel = uilabel(app.rocket);
            app.InjectiontimesLabel.HorizontalAlignment = 'right';
            app.InjectiontimesLabel.FontSize = 8;
            app.InjectiontimesLabel.Interpreter = 'tex';
            app.InjectiontimesLabel.Position = [111 48 72 22];
            app.InjectiontimesLabel.Text = 'Injection time/s';

            % Create thrust_time
            app.thrust_time = uieditfield(app.rocket, 'numeric');
            app.thrust_time.ValueChangedFcn = createCallbackFcn(app, @thrust_timeValueChanged, true);
            app.thrust_time.FontSize = 8;
            app.thrust_time.Position = [212 46 36 27];
            app.thrust_time.Value = 1.6;

            % Create Label_11
            app.Label_11 = uilabel(app.rocket);
            app.Label_11.HorizontalAlignment = 'right';
            app.Label_11.FontSize = 8;
            app.Label_11.Interpreter = 'tex';
            app.Label_11.Position = [263 61 114 22];
            app.Label_11.Text = 'Fuel mass at ignition';

            % Create fuel_weight_begin
            app.fuel_weight_begin = uieditfield(app.rocket, 'numeric');
            app.fuel_weight_begin.ValueChangedFcn = createCallbackFcn(app, @fuel_weight_beginValueChanged, true);
            app.fuel_weight_begin.FontSize = 8;
            app.fuel_weight_begin.Position = [399 60 38 25];
            app.fuel_weight_begin.Value = 25.8;

            % Create NsecLabel
            app.NsecLabel = uilabel(app.rocket);
            app.NsecLabel.HorizontalAlignment = 'right';
            app.NsecLabel.FontSize = 8;
            app.NsecLabel.Interpreter = 'tex';
            app.NsecLabel.Position = [93 17 90 22];
            app.NsecLabel.Text = 'Total Inpulse/N*sec';

            % Create total_impulse
            app.total_impulse = uieditfield(app.rocket, 'numeric');
            app.total_impulse.ValueChangedFcn = createCallbackFcn(app, @total_impulseValueChanged, true);
            app.total_impulse.FontSize = 8;
            app.total_impulse.Position = [212 14 36 28];
            app.total_impulse.Value = 10;

            % Create kgLabel
            app.kgLabel = uilabel(app.rocket);
            app.kgLabel.HorizontalAlignment = 'right';
            app.kgLabel.FontSize = 8;
            app.kgLabel.Interpreter = 'tex';
            app.kgLabel.Position = [266 31 122 22];
            app.kgLabel.Text = ' Fuel mass at end of burn ';

            % Create fuel_weight_end
            app.fuel_weight_end = uieditfield(app.rocket, 'numeric');
            app.fuel_weight_end.ValueChangedFcn = createCallbackFcn(app, @fuel_weight_endValueChanged2, true);
            app.fuel_weight_end.FontSize = 8;
            app.fuel_weight_end.Position = [399 33 38 22];
            app.fuel_weight_end.Value = 12.48;

            % Create Label_12
            app.Label_12 = uilabel(app.rocket);
            app.Label_12.HorizontalAlignment = 'right';
            app.Label_12.FontSize = 8;
            app.Label_12.Interpreter = 'tex';
            app.Label_12.Position = [254 118 157 22];
            app.Label_12.Text = 'Cross-sectional area of fuselage/m^{2}';

            % Create cross_sectional_side_area
            app.cross_sectional_side_area = uieditfield(app.rocket, 'numeric');
            app.cross_sectional_side_area.ValueChangedFcn = createCallbackFcn(app, @cross_sectional_side_areaValueChanged, true);
            app.cross_sectional_side_area.FontSize = 8;
            app.cross_sectional_side_area.Position = [418 113 41 29];
            app.cross_sectional_side_area.Value = 0.001;

            % Create WeightunitsettingLabel
            app.WeightunitsettingLabel = uilabel(app.rocket);
            app.WeightunitsettingLabel.FontSize = 8;
            app.WeightunitsettingLabel.Interpreter = 'tex';
            app.WeightunitsettingLabel.Position = [305 10 88 22];
            app.WeightunitsettingLabel.Text = 'Weight unit setting';

            % Create figures
            app.figures = uitabgroup(app.UIFigure);
            app.figures.Position = [554 297 406 322];

            % Create ThefinaldestinationTab
            app.ThefinaldestinationTab = uitab(app.figures);
            app.ThefinaldestinationTab.Title = 'The final destination';

            % Create final_destination
            app.final_destination = uiaxes(app.ThefinaldestinationTab);
            title(app.final_destination, 'The final destination')
            xlabel(app.final_destination, 'West')
            ylabel(app.final_destination, 'North')
            zlabel(app.final_destination, 'Z')
            app.final_destination.XGrid = 'on';
            app.final_destination.YGrid = 'on';
            app.final_destination.Position = [5 3 401 291];

            % Create OrbitTab
            app.OrbitTab = uitab(app.figures);
            app.OrbitTab.Title = 'Orbit';

            % Create orbit
            app.orbit = uiaxes(app.OrbitTab);
            title(app.orbit, 'Orbit')
            xlabel(app.orbit, 'West')
            ylabel(app.orbit, 'North')
            zlabel(app.orbit, 'Z')
            app.orbit.XAxisLocation = 'origin';
            app.orbit.XTick = [0 0.2 0.4 0.6 0.8 1];
            app.orbit.XTickLabel = {''; '0.2'; '0.4'; ''; '0.8'; '1'};
            app.orbit.XGrid = 'on';
            app.orbit.YGrid = 'on';
            app.orbit.ZGrid = 'on';
            app.orbit.ButtonDownFcn = createCallbackFcn(app, @ButtonValueChanged, true);
            app.orbit.Position = [8 13 389 281];

            % Create vtgraphTab
            app.vtgraphTab = uitab(app.figures);
            app.vtgraphTab.Title = 'v-t graph';

            % Create v_t_graph
            app.v_t_graph = uiaxes(app.vtgraphTab);
            title(app.v_t_graph, 'v-t graph')
            xlabel(app.v_t_graph, 'TIme/s')
            ylabel(app.v_t_graph, 'Velocity/ms^{-1}')
            zlabel(app.v_t_graph, 'Z')
            app.v_t_graph.Position = [8 3 398 283];

            % Create htgraphTab
            app.htgraphTab = uitab(app.figures);
            app.htgraphTab.Title = 'h-t graph';

            % Create h_t_graph
            app.h_t_graph = uiaxes(app.htgraphTab);
            title(app.h_t_graph, 'h-t graph')
            xlabel(app.h_t_graph, 'TIme/s')
            ylabel(app.h_t_graph, 'Height/m')
            zlabel(app.h_t_graph, 'Z')
            app.h_t_graph.Position = [2 3 403 291];

            % Create vhtTab
            app.vhtTab = uitab(app.figures);
            app.vhtTab.Title = 'v-h-t';

            % Create v_h_t_graph
            app.v_h_t_graph = uiaxes(app.vhtTab);
            title(app.v_h_t_graph, 'v-h-t')
            xlabel(app.v_h_t_graph, 'X')
            ylabel(app.v_h_t_graph, 'Y')
            zlabel(app.v_h_t_graph, 'Z')
            app.v_h_t_graph.Position = [5 3 389 289];

            % Create Panel_parachute
            app.Panel_parachute = uipanel(app.UIFigure);
            app.Panel_parachute.Title = 'Parachute setting';
            app.Panel_parachute.Position = [15 185 304 256];

            % Create parachute_on_off
            app.parachute_on_off = uidropdown(app.Panel_parachute);
            app.parachute_on_off.Items = {'Parachute on', 'Parachute off'};
            app.parachute_on_off.ValueChangedFcn = createCallbackFcn(app, @parachute_on_offValueChanged, true);
            app.parachute_on_off.Position = [157 189 136 35];
            app.parachute_on_off.Value = 'Parachute on';

            % Create Label_5
            app.Label_5 = uilabel(app.Panel_parachute);
            app.Label_5.HorizontalAlignment = 'right';
            app.Label_5.FontSize = 8;
            app.Label_5.Interpreter = 'tex';
            app.Label_5.Position = [17 138 199 22];
            app.Label_5.Text = 'Time to parachute deployment/s';

            % Create parachute_time
            app.parachute_time = uieditfield(app.Panel_parachute, 'numeric');
            app.parachute_time.ValueChangedFcn = createCallbackFcn(app, @parachute_timeValueChanged, true);
            app.parachute_time.FontSize = 8;
            app.parachute_time.Position = [226 138 39 27];
            app.parachute_time.Value = 3;

            % Create Label_6
            app.Label_6 = uilabel(app.Panel_parachute);
            app.Label_6.HorizontalAlignment = 'right';
            app.Label_6.FontSize = 8;
            app.Label_6.Interpreter = 'tex';
            app.Label_6.Position = [22 100 180 22];
            app.Label_6.Text = 'Projected area of parachute/m^{2}';

            % Create parachute_S
            app.parachute_S = uieditfield(app.Panel_parachute, 'numeric');
            app.parachute_S.ValueChangedFcn = createCallbackFcn(app, @parachute_SValueChanged, true);
            app.parachute_S.FontSize = 8;
            app.parachute_S.Position = [213 97 52 28];
            app.parachute_S.Value = 0.1;

            % Create Label_7
            app.Label_7 = uilabel(app.Panel_parachute);
            app.Label_7.HorizontalAlignment = 'right';
            app.Label_7.FontSize = 8;
            app.Label_7.Interpreter = 'tex';
            app.Label_7.Position = [20 56 186 31];
            app.Label_7.Text = 'Resistance coefficient of parachute(Cd)';

            % Create drag_coefficient_parachute
            app.drag_coefficient_parachute = uieditfield(app.Panel_parachute, 'numeric');
            app.drag_coefficient_parachute.ValueChangedFcn = createCallbackFcn(app, @drag_coefficient_parachuteValueChanged, true);
            app.drag_coefficient_parachute.FontSize = 8;
            app.drag_coefficient_parachute.Position = [213 57 52 30];
            app.drag_coefficient_parachute.Value = 3;

            % Create Label_8
            app.Label_8 = uilabel(app.Panel_parachute);
            app.Label_8.HorizontalAlignment = 'right';
            app.Label_8.FontSize = 8;
            app.Label_8.Interpreter = 'tex';
            app.Label_8.Position = [32 14 172 31];
            app.Label_8.Text = 'Parachute canopy load coefficient(Cxo)';

            % Create C_xo
            app.C_xo = uieditfield(app.Panel_parachute, 'numeric');
            app.C_xo.ValueChangedFcn = createCallbackFcn(app, @C_xoValueChanged, true);
            app.C_xo.FontSize = 8;
            app.C_xo.Position = [213 14 52 30];
            app.C_xo.Value = 2;

            % Create fundamental_details
            app.fundamental_details = uipanel(app.UIFigure);
            app.fundamental_details.Title = 'Basic setting';
            app.fundamental_details.Position = [318 185 222 256];

            % Create Label_3
            app.Label_3 = uilabel(app.fundamental_details);
            app.Label_3.HorizontalAlignment = 'right';
            app.Label_3.FontSize = 8;
            app.Label_3.Interpreter = 'tex';
            app.Label_3.Position = [25 196 101 22];
            app.Label_3.Text = 'Air Density/kgm^{-3}';

            % Create rho
            app.rho = uieditfield(app.fundamental_details, 'numeric');
            app.rho.ValueChangedFcn = createCallbackFcn(app, @rhoValueChanged, true);
            app.rho.FontSize = 8;
            app.rho.Position = [137 190 65 33];
            app.rho.Value = 1.225;

            % Create GravitationalaccelerationLabel
            app.GravitationalaccelerationLabel = uilabel(app.fundamental_details);
            app.GravitationalaccelerationLabel.HorizontalAlignment = 'right';
            app.GravitationalaccelerationLabel.FontSize = 8;
            app.GravitationalaccelerationLabel.Interpreter = 'tex';
            app.GravitationalaccelerationLabel.Position = [25 155 108 22];
            app.GravitationalaccelerationLabel.Text = 'Gravitational acceleration';

            % Create gravitational_acceleration
            app.gravitational_acceleration = uieditfield(app.fundamental_details, 'numeric');
            app.gravitational_acceleration.ValueChangedFcn = createCallbackFcn(app, @gravitational_accelerationValueChanged, true);
            app.gravitational_acceleration.FontSize = 8;
            app.gravitational_acceleration.Position = [138 149 64 34];
            app.gravitational_acceleration.Value = 9.8;

            % Create Airvelocityms1Label
            app.Airvelocityms1Label = uilabel(app.fundamental_details);
            app.Airvelocityms1Label.HorizontalAlignment = 'right';
            app.Airvelocityms1Label.FontSize = 8;
            app.Airvelocityms1Label.Interpreter = 'tex';
            app.Airvelocityms1Label.Position = [24 112 89 22];
            app.Airvelocityms1Label.Text = 'Wind velocity/ms^{-1}';

            % Create wind_speed
            app.wind_speed = uieditfield(app.fundamental_details, 'numeric');
            app.wind_speed.ValueChangedFcn = createCallbackFcn(app, @wind_speedValueChanged, true);
            app.wind_speed.FontSize = 8;
            app.wind_speed.Position = [137 106 65 34];
            app.wind_speed.Value = 1.5;

            % Create WinddirectiondegLabel
            app.WinddirectiondegLabel = uilabel(app.fundamental_details);
            app.WinddirectiondegLabel.HorizontalAlignment = 'right';
            app.WinddirectiondegLabel.FontSize = 8;
            app.WinddirectiondegLabel.Interpreter = 'tex';
            app.WinddirectiondegLabel.Position = [26 71 87 22];
            app.WinddirectiondegLabel.Text = 'Wind direction/deg';

            % Create wind_degree
            app.wind_degree = uieditfield(app.fundamental_details, 'numeric');
            app.wind_degree.ValueChangedFcn = createCallbackFcn(app, @wind_degreeValueChanged, true);
            app.wind_degree.FontSize = 8;
            app.wind_degree.Position = [137 66 65 32];

            % Create Label_10
            app.Label_10 = uilabel(app.fundamental_details);
            app.Label_10.HorizontalAlignment = 'right';
            app.Label_10.FontSize = 8;
            app.Label_10.Interpreter = 'tex';
            app.Label_10.Position = [60 27 61 22];
            app.Label_10.Text = 'Wind impact';

            % Create wind_effectness
            app.wind_effectness = uieditfield(app.fundamental_details, 'numeric');
            app.wind_effectness.ValueChangedFcn = createCallbackFcn(app, @wind_effectnessValueChanged, true);
            app.wind_effectness.FontSize = 8;
            app.wind_effectness.Position = [138 24 63 31];
            app.wind_effectness.Value = 0.1;

            % Create simulation_details
            app.simulation_details = uipanel(app.UIFigure);
            app.simulation_details.Title = 'Simulation setting';
            app.simulation_details.Position = [380 8 250 134];

            % Create SimulationtimesLabel
            app.SimulationtimesLabel = uilabel(app.simulation_details);
            app.SimulationtimesLabel.HorizontalAlignment = 'right';
            app.SimulationtimesLabel.FontSize = 8;
            app.SimulationtimesLabel.Interpreter = 'tex';
            app.SimulationtimesLabel.Position = [17 63 157 22];
            app.SimulationtimesLabel.Text = 'Simulation time/s';

            % Create simulation_time
            app.simulation_time = uieditfield(app.simulation_details, 'numeric');
            app.simulation_time.ValueChangedFcn = createCallbackFcn(app, @simulation_timeValueChanged, true);
            app.simulation_time.FontSize = 8;
            app.simulation_time.Position = [180 57 46 34];
            app.simulation_time.Value = 100;

            % Create Label
            app.Label = uilabel(app.simulation_details);
            app.Label.HorizontalAlignment = 'right';
            app.Label.FontSize = 8;
            app.Label.Interpreter = 'tex';
            app.Label.Position = [47 17 92 22];
            app.Label.Text = 'Number of attempts';

            % Create num_simulations
            app.num_simulations = uieditfield(app.simulation_details, 'numeric');
            app.num_simulations.ValueChangedFcn = createCallbackFcn(app, @num_simulationsValueChanged, true);
            app.num_simulations.FontSize = 8;
            app.num_simulations.Position = [153 11 73 34];
            app.num_simulations.Value = 1;

            % Create InitializeallgraphsButton
            app.InitializeallgraphsButton = uibutton(app.UIFigure, 'push');
            app.InitializeallgraphsButton.ButtonPushedFcn = createCallbackFcn(app, @InitializeallgraphsButtonPushed, true);
            app.InitializeallgraphsButton.Position = [831 234 129 49];
            app.InitializeallgraphsButton.Text = 'Initialize all graphs';

            % Create Panel_launch_deg
            app.Panel_launch_deg = uipanel(app.UIFigure);
            app.Panel_launch_deg.Title = 'Launch direction setting';
            app.Panel_launch_deg.Position = [15 7 365 175];

            % Create ZdegLabel
            app.ZdegLabel = uilabel(app.Panel_launch_deg);
            app.ZdegLabel.HorizontalAlignment = 'right';
            app.ZdegLabel.FontSize = 8;
            app.ZdegLabel.Interpreter = 'tex';
            app.ZdegLabel.Position = [11 99 123 22];
            app.ZdegLabel.Text = 'Launch direction \theta/deg';

            % Create theta
            app.theta = uieditfield(app.Panel_launch_deg, 'numeric');
            app.theta.ValueChangedFcn = createCallbackFcn(app, @thetaValueChanged, true);
            app.theta.FontSize = 8;
            app.theta.Position = [140 93 73 34];

            % Create Image
            app.Image = uiimage(app.Panel_launch_deg);
            app.Image.ScaleMethod = 'stretch';
            app.Image.Position = [225 17 123 110];
            app.Image.ImageSource = fullfile(pathToMLAPP, '1235_0_rthetaphi1.png');

            % Create LaunchdirectionphidegLabel
            app.LaunchdirectionphidegLabel = uilabel(app.Panel_launch_deg);
            app.LaunchdirectionphidegLabel.HorizontalAlignment = 'right';
            app.LaunchdirectionphidegLabel.FontSize = 8;
            app.LaunchdirectionphidegLabel.Interpreter = 'tex';
            app.LaunchdirectionphidegLabel.Position = [10 58 125 22];
            app.LaunchdirectionphidegLabel.Text = 'Launch direction \phi/deg';

            % Create phi
            app.phi = uieditfield(app.Panel_launch_deg, 'numeric');
            app.phi.ValueChangedFcn = createCallbackFcn(app, @phiValueChanged, true);
            app.phi.FontSize = 8;
            app.phi.Position = [141 52 73 34];

            % Create start_button
            app.start_button = uibutton(app.UIFigure, 'state');
            app.start_button.ValueChangedFcn = createCallbackFcn(app, @ButtonValueChanged, true);
            app.start_button.Text = 'Simulation start';
            app.start_button.Position = [748 181 212 48];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = rocket_orbit_simulator_test_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end