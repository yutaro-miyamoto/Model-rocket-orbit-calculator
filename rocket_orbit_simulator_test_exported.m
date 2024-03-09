classdef rocket_orbit_simulator_test_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        Button_2                    matlab.ui.control.Button
        Panel_launch_deg            matlab.ui.container.Panel
        Image                       matlab.ui.control.Image
        phi                         matlab.ui.control.NumericEditField
        phidegEditFieldLabel        matlab.ui.control.Label
        theta                       matlab.ui.control.NumericEditField
        ZdegLabel                   matlab.ui.control.Label
        rocket                      matlab.ui.container.Panel
        checking_kg_or_g            matlab.ui.control.DropDown
        fuel_weight_begin           matlab.ui.control.NumericEditField
        Label_11                    matlab.ui.control.Label
        fuel_weight_end             matlab.ui.control.NumericEditField
        kgLabel                     matlab.ui.control.Label
        S_vertical                  matlab.ui.control.NumericEditField
        Label_4                     matlab.ui.control.Label
        drag_coefficient            matlab.ui.control.NumericEditField
        Label_2                     matlab.ui.control.Label
        thrust_time                 matlab.ui.control.NumericEditField
        sEditFieldLabel_2           matlab.ui.control.Label
        total_impulse               matlab.ui.control.NumericEditField
        NsecLabel                   matlab.ui.control.Label
        constant_weight             matlab.ui.control.NumericEditField
        gLabel                      matlab.ui.control.Label
        Panel_parachute             matlab.ui.container.Panel
        parachute_on_off            matlab.ui.control.DropDown
        C_xo                        matlab.ui.control.NumericEditField
        Label_8                     matlab.ui.control.Label
        drag_coefficient_parachute  matlab.ui.control.NumericEditField
        Label_7                     matlab.ui.control.Label
        parachute_S                 matlab.ui.control.NumericEditField
        Label_6                     matlab.ui.control.Label
        parachute_time              matlab.ui.control.NumericEditField
        Label_5                     matlab.ui.control.Label
        simulation_details          matlab.ui.container.Panel
        num_simulations             matlab.ui.control.NumericEditField
        Label                       matlab.ui.control.Label
        simulation_time             matlab.ui.control.NumericEditField
        sEditFieldLabel             matlab.ui.control.Label
        fundamental_details         matlab.ui.container.Panel
        wind_effectness             matlab.ui.control.NumericEditField
        Label_10                    matlab.ui.control.Label
        wind_degree                 matlab.ui.control.NumericEditField
        degEditFieldLabel           matlab.ui.control.Label
        rho                         matlab.ui.control.NumericEditField
        Label_3                     matlab.ui.control.Label
        wind_speed                  matlab.ui.control.NumericEditField
        ms1EditFieldLabel           matlab.ui.control.Label
        gravitational_acceleration  matlab.ui.control.NumericEditField
        EditFieldLabel              matlab.ui.control.Label
        figures                     matlab.ui.container.TabGroup
        Tab                         matlab.ui.container.Tab
        final_destination           matlab.ui.control.UIAxes
        Tab_2                       matlab.ui.container.Tab
        orbit                       matlab.ui.control.UIAxes
        Tab_3                       matlab.ui.container.Tab
        v_t_graph                   matlab.ui.control.UIAxes
        Tab_4                       matlab.ui.container.Tab
        h_t_graph                   matlab.ui.control.UIAxes
        Tab_5                       matlab.ui.container.Tab
        v_h_t_graph                 matlab.ui.control.UIAxes
        start_button                matlab.ui.control.StateButton
    end

   

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % 外部関数があるディレクトリを追加
        addpath('/functions');
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

        % Value changed function: S_vertical
        function S_verticalValueChanged(app, event)
            % 機体断面積
            value = app.S_vertical.Value;            
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

        % Button pushed function: Button_2
        function Button_2Pushed(app, event)
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
                tmax = 100.0;

                %kg or g の決定
                if strcmpi(app.checking_kg_or_g.Value, 'kg')
                    times_const = 1;
                elseif strcmpi(app.checking_kg_or_g.Value, 'g')
                    times_const = 0.001;
                end
                %パラシュートのon/offを切り替え
                if strcmpi(app.parachute_on_off.Value,'パラシュート無し')
                    app.parachute_time.Value = tmax;
                end

                % シミュレーションの開始
                num_simulations = app.num_simulations.Value;

                % 最終到着地点を保存する配列
                landing_points = zeros(2, num_simulations);
                landing_points_1 = zeros(2, num_simulations);
                for a = 1:num_simulations
                    % 初期位置
                    r0 = [0;0;0];
                    r_1_0 = r0;

                    %　総消費燃料質量(消費される推進剤の全質量)
                    m_fuel = app.fuel_weight_begin.Value - app.fuel_weight_end.Value;
                    m_fuel = m_fuel * times_const;%kgへの単位変換

                    % 全質量/g
                    m = app.constant_weight.Value + app.fuel_weight_begin.Value;
                    m = m * times_const;%kgへの単位変換

                    % 時間ステップ/s
                    dt = 0.001;

                    %(dt)s当たりの消費燃料質量
                    m_fuel_using = m_fuel / app.thrust_time.Value * dt;


                    %推力/N
                    F_r0 = [0;0;app.total_impulse.Value / app.thrust_time.Value];

                    %平均比推力Isp
                    Isp = app.total_impulse.Value / m_fuel / app.gravitational_acceleration.Value;

                    % 初期速度を打ち上げ角度に基づいて計算
                    %v0 = [0;0;0];
                    v_1_0 = [0;0;0];
                    %v = v0;
                    v_1 = v_1_0;
                    %sv=v0.*v0;% ベクトルの二乗
                    %dp=sum(sv);% 二乗和
                    %mag_v=sqrt(dp);% 絶対値

                    heading_vector = [sind(app.theta.Value)*cosd(app.phi.Value);sind(app.theta.Value)*sind(app.phi.Value);cosd(app.theta.Value)];
                    F_rv = F_r0(3) * heading_vector;

                    %平均排気速度v
                    Dv = app.total_impulse.Value / m_fuel * heading_vector;

                    % 軌道データを保存する配列
                    %r_data = zeros(3, tmax / dt);
                    r_1_data = zeros(3, tmax / dt);

                    % 初期位置を保存
                    %r_data(:,1) = r0;
                    r_1_data(:,1) = r0;

                    % 時間データを保存する配列
                    t_data = zeros(1, tmax / dt);

                    % 速度データを保存する配列
                    %v_data = zeros(1, tmax / dt);
                    v_1_data = zeros(1, tmax / dt);

                    % ループを実行
                    % 風速度
                    %wind_speed = [app.wind_speed.Value*cosd(app.wind_degree.Value) + app.wind_effectness.Value * randn;...
                    %    app.wind_speed.Value*sind(app.wind_degree.Value)+ app.wind_effectness.Value *randn;0];
                    wind_speed = [0;0;0];%テスト用

                    for t = 1:tmax/dt   %注意！！！この中の計算はすべてdt秒あたり、で考えないと値がおかしくなる。

                        %空気抵抗の計算
                        %Air_resistance = app.drag_coefficient.Value * 0.5 * app.rho.Value * sqrt(sum(v.*v)) * v * app.S_vertical.Value;
                        Air_resistance_1 = app.drag_coefficient.Value * 0.5 * app.rho.Value * sqrt(sum(v_1.^2)) * v_1 * app.S_vertical.Value;


                        if t <= app.thrust_time.Value / dt %%ロケットついてるとき
                            % 力
                            %F = [0;0;-m*app.gravitational_acceleration.Value] + F_rv - Air_resistance;
                            F_1 = calculate_force(m, app.gravitational_acceleration.Value, Air_resistance_1, F_rv);
                            %F_1 = [0;0;-m*app.gravitational_acceleration.Value] + F_rv - Air_resistance_1;

                            % 速度を更新
                            %v = v0 + F / m * dt;
                            v_1 = velocity(v_1_0,Dv,m_fuel_using,m,F_1,dt);
                            %v_1 =  v_1_0 + Dv * m_fuel_using / m +  F_1 / m * dt;

                            

                            % 位置を更新
                            %r = r0 + (v + wind_speed) * dt;
                            r_1 = location(r_1_0, v_1, wind_speed, dt);
                            %r_1 = r_1_0 + (v_1 + wind_speed) * dt;

                            % 結果を保存
                            %r_data(:,t) = r;
                            r_1_data(:,t) = r_1;

                            % 質量の更新(燃料の消費)
                            m1 = m;
                            m = m1 - m_fuel_using * dt;%dtあたりの燃料消費量を引く

                        elseif  t > app.thrust_time.Value / dt && t <= (app.thrust_time.Value + app.parachute_time.Value) / dt
                            %%ロケット推進消えたE-

                            % 力
                            %F = [0;0;-m * app.gravitational_acceleration.Value] - Air_resistance;
                            F_1 = [0;0;-m * app.gravitational_acceleration.Value] - Air_resistance_1;

                            % 速度を更新
                            %v = v0 + F / m * dt;
                            v_1 = v_1_0 + F_1 / m * dt;


                            % 位置を更新
                            %r = r0 + (v + wind_speed) * dt;
                            r_1 = r_1_0 + (v_1 + wind_speed) * dt;

                            % 結果を保存
                            %r_data(:,t) = r;
                            r_1_data(:,t) = r_1;

                        else%パラシュート展開。抗力係数と機体の投影面積が変わる。
                            %空気抵抗の更新
                            %https://nociws.github.io/parachute/
                            %を参照
                            %Air_resistance = app.drag_coefficient_parachute.Value * 0.5 * app.rho.Value * sqrt(sum(v.^2))...
                            %    * v * app.parachute_S.Value * app.C_xo.Value;
                            Air_resistance_1 = app.drag_coefficient_parachute.Value * 0.5 * app.rho.Value * sqrt(sum(v_1.^2))...
                                * v_1 * app.parachute_S.Value * app.C_xo.Value;

                            % 力
                            %F = [0;0;-m*app.gravitational_acceleration.Value] - Air_resistance;
                            F_1 = [0;0;-m * app.gravitational_acceleration.Value] - Air_resistance_1;

                            % 速度を更新
                            %v = v0 + F / m * dt;
                            v_1 = v_1_0 + F_1 / m * dt;

                            % 位置を更新
                            %r = r0 + (v + wind_speed) * dt;
                            r_1 = r_1_0 + (v_1 + wind_speed) * dt;

                            % 結果を保存
                            %r_data(:,t) = r;
                            r_1_data(:,t) = r_1;

                        end
                        % rのz座標がゼロになった場合、ループから抜ける
                        %if r(3)<0
                        % z座標が無視されることを保証するために、rの3番目の要素を0にする
                        %r(3) = 0;

                        % 結果を保存
                        %landing_points(:,a) = r(1:2); % x座標とy座標のみを保存する
                        %break;
                        if r_1(3)<0                            
                            r_1(3) = 0;
                            landing_points_1(:,a) = r_1(1:2);
                            break
                        end
                        %end

                        % 次の時間ステップの準備
                        % 速度ベクトルを用いたヘディングベクトルの生成
                        %mag_v = sqrt(sum(v.^2));% 絶対値
                        mag_1_v = sqrt(sum(v_1.^2));

                        %v_data(t) = mag_v;
                        v_1_data(t) = mag_1_v;
                        t_data(t) = t*dt;

                        %heading_vector = v / mag_v; %これ上のif文にいれないで！！
                        heading_vector_1 = v_1 / mag_1_v;
                        F_rv = F_r0(3) * heading_vector;%動かなくなるよ～

                        Dv = app.total_impulse.Value / m_fuel * heading_vector_1;

                        %r0 = r;
                        r_1_0= r_1;

                        %v0 = v;
                        v_1_0 = v_1;

                    end

                    %plot3(app.orbit, r_data(1,:), r_data(2,:), r_data(3,:));
                    %hold(app.orbit, 'on');
                    plot3(app.orbit, r_1_data(1,:), r_1_data(2,:), r_1_data(3,:));
                    hold(app.orbit, 'on'); % Axes コンポーネントに hold on を設定

                    %plot(app.v_t_graph, t_data(1,:), v_data(1,:));
                    %hold(app.v_t_graph, 'on');
                    plot(app.v_t_graph, t_data(1,:), v_1_data(1,:));
                    hold(app.v_t_graph, 'on');

                    % plot(app.h_t_graph, t_data(1,:), r_data(3,:));
                    % hold(app.h_t_graph, 'on');
                    plot(app.h_t_graph, t_data(1,:), r_1_data(3,:));
                    hold(app.h_t_graph, 'on');

                    yyaxis(app.v_h_t_graph, 'left');
                    % plot(app.v_h_t_graph, t_data(1,:), v_data(1,:));
                    % hold(app.v_h_t_graph, 'on');
                    plot(app.v_h_t_graph, t_data(1,:), v_1_data(1,:));
                    hold(app.v_h_t_graph, 'on');
                    % if max(v_1_data(1,:))>max(v_data(1,:))
                    %     ylim(app.v_h_t_graph, [0, max(v_1_data(1,:))]); % 左の軸の範囲を設定
                    % else
                    %     ylim(app.v_h_t_graph, [0, max(v_data(1,:))]);
                    % end
                    ylim(app.v_h_t_graph, [0, max(v_1_data(1,:))]);

                    yyaxis(app.v_h_t_graph, 'right');
                    % plot(app.v_h_t_graph, t_data(1,:), r_data(3,:));
                    % hold(app.v_h_t_graph, "on");
                    plot(app.v_h_t_graph, t_data(1,:), r_1_data(3,:));
                    hold(app.v_h_t_graph, "on");
                    % if max(r_data(3,:)) > max(r_1_data(3,:))
                    %     ylim(app.v_h_t_graph, [0, max(r_data(3,:))]); % 右の軸の範囲を設定
                    % else
                    %     ylim(app.v_h_t_graph, [0, max(r_1_data(3,:))]);
                    % end
                    ylim(app.v_h_t_graph, [0, max(r_1_data(3,:))]);

                    ylabel(app.v_h_t_graph, '速度');
                    ylabel(app.v_h_t_graph, '高度');

                end
                % シミュレーションの結果をプロットするコード
                %scatter(app.final_destination, landing_points(1,:), landing_points(2,:), 'magenta');
                scatter(app.final_destination, landing_points_1(1,:), landing_points_1(2,:), 'yellow');
                hold(app.orbit, 'off'); % Axes コンポーネントに hold off を設定
                hold(app.final_destination, 'off'); % Axes コンポーネントに hold off を設定
                hold(app.v_t_graph, 'off'); % Axes コンポーネントに hold off を設定
                hold(app.h_t_graph, 'off'); % Axes コンポーネントに hold off を設定
                hold(app.v_h_t_graph, 'off'); % Axes コンポーネントに hold off を設定

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
            app.UIFigure.Position = [100 100 1086 685];
            app.UIFigure.Name = 'MATLAB App';

            % Create start_button
            app.start_button = uibutton(app.UIFigure, 'state');
            app.start_button.ValueChangedFcn = createCallbackFcn(app, @ButtonValueChanged, true);
            app.start_button.Text = 'シュミレーションスタート';
            app.start_button.Position = [677 21 212 48];

            % Create figures
            app.figures = uitabgroup(app.UIFigure);
            app.figures.Position = [556 284 523 373];

            % Create Tab
            app.Tab = uitab(app.figures);
            app.Tab.Title = '最終到着地点';

            % Create final_destination
            app.final_destination = uiaxes(app.Tab);
            title(app.final_destination, '最終到着地点(xy平面)')
            xlabel(app.final_destination, '東')
            ylabel(app.final_destination, '北')
            zlabel(app.final_destination, 'Z')
            app.final_destination.XGrid = 'on';
            app.final_destination.YGrid = 'on';
            app.final_destination.Position = [5 7 513 337];

            % Create Tab_2
            app.Tab_2 = uitab(app.figures);
            app.Tab_2.Title = '軌道';

            % Create orbit
            app.orbit = uiaxes(app.Tab_2);
            title(app.orbit, '軌道')
            xlabel(app.orbit, '東')
            ylabel(app.orbit, '北')
            zlabel(app.orbit, 'Z')
            app.orbit.XAxisLocation = 'origin';
            app.orbit.XTick = [0 0.2 0.4 0.6 0.8 1];
            app.orbit.XTickLabel = {''; '0.2'; '0.4'; ''; '0.8'; '1'};
            app.orbit.XGrid = 'on';
            app.orbit.YGrid = 'on';
            app.orbit.ZGrid = 'on';
            app.orbit.ButtonDownFcn = createCallbackFcn(app, @ButtonValueChanged, true);
            app.orbit.Position = [18 7 486 328];

            % Create Tab_3
            app.Tab_3 = uitab(app.figures);
            app.Tab_3.Title = '速度/時間';

            % Create v_t_graph
            app.v_t_graph = uiaxes(app.Tab_3);
            title(app.v_t_graph, 'v-t graph')
            xlabel(app.v_t_graph, '時間')
            ylabel(app.v_t_graph, '速度')
            zlabel(app.v_t_graph, 'Z')
            app.v_t_graph.Position = [8 7 510 330];

            % Create Tab_4
            app.Tab_4 = uitab(app.figures);
            app.Tab_4.Title = '高度/時間';

            % Create h_t_graph
            app.h_t_graph = uiaxes(app.Tab_4);
            title(app.h_t_graph, 'h-t graph')
            xlabel(app.h_t_graph, '時間')
            ylabel(app.h_t_graph, '高度')
            zlabel(app.h_t_graph, 'Z')
            app.h_t_graph.Position = [3 7 515 336];

            % Create Tab_5
            app.Tab_5 = uitab(app.figures);
            app.Tab_5.Title = '高度/速度/時間';

            % Create v_h_t_graph
            app.v_h_t_graph = uiaxes(app.Tab_5);
            title(app.v_h_t_graph, 'v-h-t')
            xlabel(app.v_h_t_graph, 'X')
            ylabel(app.v_h_t_graph, 'Y')
            zlabel(app.v_h_t_graph, 'Z')
            app.v_h_t_graph.Position = [5 12 513 331];

            % Create fundamental_details
            app.fundamental_details = uipanel(app.UIFigure);
            app.fundamental_details.Title = '基礎物理量';
            app.fundamental_details.Position = [313 142 230 313];

            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.fundamental_details);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Interpreter = 'tex';
            app.EditFieldLabel.Position = [35 212 78 22];
            app.EditFieldLabel.Text = '重力加速度';

            % Create gravitational_acceleration
            app.gravitational_acceleration = uieditfield(app.fundamental_details, 'numeric');
            app.gravitational_acceleration.ValueChangedFcn = createCallbackFcn(app, @gravitational_accelerationValueChanged, true);
            app.gravitational_acceleration.Position = [137 206 65 34];
            app.gravitational_acceleration.Value = 9.8;

            % Create ms1EditFieldLabel
            app.ms1EditFieldLabel = uilabel(app.fundamental_details);
            app.ms1EditFieldLabel.HorizontalAlignment = 'right';
            app.ms1EditFieldLabel.Interpreter = 'tex';
            app.ms1EditFieldLabel.Position = [35 169 78 22];
            app.ms1EditFieldLabel.Text = '風速度/ms^{-1}';

            % Create wind_speed
            app.wind_speed = uieditfield(app.fundamental_details, 'numeric');
            app.wind_speed.ValueChangedFcn = createCallbackFcn(app, @wind_speedValueChanged, true);
            app.wind_speed.Position = [137 163 65 34];

            % Create Label_3
            app.Label_3 = uilabel(app.fundamental_details);
            app.Label_3.HorizontalAlignment = 'right';
            app.Label_3.Interpreter = 'tex';
            app.Label_3.Position = [25 253 101 22];
            app.Label_3.Text = '空気密度/kgm^{-3}';

            % Create rho
            app.rho = uieditfield(app.fundamental_details, 'numeric');
            app.rho.ValueChangedFcn = createCallbackFcn(app, @rhoValueChanged, true);
            app.rho.Position = [137 247 65 33];
            app.rho.Value = 1.225;

            % Create degEditFieldLabel
            app.degEditFieldLabel = uilabel(app.fundamental_details);
            app.degEditFieldLabel.HorizontalAlignment = 'right';
            app.degEditFieldLabel.Interpreter = 'tex';
            app.degEditFieldLabel.Position = [39 128 74 22];
            app.degEditFieldLabel.Text = '風向き/deg';

            % Create wind_degree
            app.wind_degree = uieditfield(app.fundamental_details, 'numeric');
            app.wind_degree.ValueChangedFcn = createCallbackFcn(app, @wind_degreeValueChanged, true);
            app.wind_degree.Position = [137 123 65 32];

            % Create Label_10
            app.Label_10 = uilabel(app.fundamental_details);
            app.Label_10.HorizontalAlignment = 'right';
            app.Label_10.Position = [76 53 53 22];
            app.Label_10.Text = '風影響度';

            % Create wind_effectness
            app.wind_effectness = uieditfield(app.fundamental_details, 'numeric');
            app.wind_effectness.ValueChangedFcn = createCallbackFcn(app, @wind_effectnessValueChanged, true);
            app.wind_effectness.Position = [144 44 70 40];

            % Create simulation_details
            app.simulation_details = uipanel(app.UIFigure);
            app.simulation_details.Title = 'シュミレーション設定';
            app.simulation_details.Position = [542 142 250 134];

            % Create sEditFieldLabel
            app.sEditFieldLabel = uilabel(app.simulation_details);
            app.sEditFieldLabel.HorizontalAlignment = 'right';
            app.sEditFieldLabel.Interpreter = 'tex';
            app.sEditFieldLabel.Position = [17 63 157 22];
            app.sEditFieldLabel.Text = 'シュミレーション時間/s';

            % Create simulation_time
            app.simulation_time = uieditfield(app.simulation_details, 'numeric');
            app.simulation_time.ValueChangedFcn = createCallbackFcn(app, @simulation_timeValueChanged, true);
            app.simulation_time.Position = [180 57 46 34];
            app.simulation_time.Value = 100;

            % Create Label
            app.Label = uilabel(app.simulation_details);
            app.Label.HorizontalAlignment = 'right';
            app.Label.Interpreter = 'tex';
            app.Label.Position = [75 17 64 22];
            app.Label.Text = '試行回数';

            % Create num_simulations
            app.num_simulations = uieditfield(app.simulation_details, 'numeric');
            app.num_simulations.ValueChangedFcn = createCallbackFcn(app, @num_simulationsValueChanged, true);
            app.num_simulations.Position = [153 11 73 34];
            app.num_simulations.Value = 1;

            % Create Panel_parachute
            app.Panel_parachute = uipanel(app.UIFigure);
            app.Panel_parachute.Title = 'パラシュート設定';
            app.Panel_parachute.Position = [10 199 304 256];

            % Create Label_5
            app.Label_5 = uilabel(app.Panel_parachute);
            app.Label_5.HorizontalAlignment = 'right';
            app.Label_5.Interpreter = 'tex';
            app.Label_5.Position = [17 138 199 22];
            app.Label_5.Text = 'パラシュート展開までの時間/s';

            % Create parachute_time
            app.parachute_time = uieditfield(app.Panel_parachute, 'numeric');
            app.parachute_time.ValueChangedFcn = createCallbackFcn(app, @parachute_timeValueChanged, true);
            app.parachute_time.Position = [226 138 39 27];
            app.parachute_time.Value = 3;

            % Create Label_6
            app.Label_6 = uilabel(app.Panel_parachute);
            app.Label_6.HorizontalAlignment = 'right';
            app.Label_6.Interpreter = 'tex';
            app.Label_6.Position = [22 100 180 22];
            app.Label_6.Text = 'パラシュートの投影面積/m^{2}';

            % Create parachute_S
            app.parachute_S = uieditfield(app.Panel_parachute, 'numeric');
            app.parachute_S.ValueChangedFcn = createCallbackFcn(app, @parachute_SValueChanged, true);
            app.parachute_S.Position = [213 97 52 28];
            app.parachute_S.Value = 0.1;

            % Create Label_7
            app.Label_7 = uilabel(app.Panel_parachute);
            app.Label_7.HorizontalAlignment = 'right';
            app.Label_7.Interpreter = 'tex';
            app.Label_7.Position = [20 56 186 31];
            app.Label_7.Text = 'パラシュートの抵抗係数(Cd)';

            % Create drag_coefficient_parachute
            app.drag_coefficient_parachute = uieditfield(app.Panel_parachute, 'numeric');
            app.drag_coefficient_parachute.ValueChangedFcn = createCallbackFcn(app, @drag_coefficient_parachuteValueChanged, true);
            app.drag_coefficient_parachute.Position = [213 57 52 30];
            app.drag_coefficient_parachute.Value = 3;

            % Create Label_8
            app.Label_8 = uilabel(app.Panel_parachute);
            app.Label_8.HorizontalAlignment = 'right';
            app.Label_8.Interpreter = 'tex';
            app.Label_8.Position = [78 14 126 31];
            app.Label_8.Text = '傘体荷重係数(Cxo)';

            % Create C_xo
            app.C_xo = uieditfield(app.Panel_parachute, 'numeric');
            app.C_xo.ValueChangedFcn = createCallbackFcn(app, @C_xoValueChanged, true);
            app.C_xo.Position = [213 14 52 30];
            app.C_xo.Value = 2;

            % Create parachute_on_off
            app.parachute_on_off = uidropdown(app.Panel_parachute);
            app.parachute_on_off.Items = {'パラシュートあり', 'パラシュート無し'};
            app.parachute_on_off.ValueChangedFcn = createCallbackFcn(app, @parachute_on_offValueChanged, true);
            app.parachute_on_off.Position = [157 189 136 35];
            app.parachute_on_off.Value = 'パラシュート無し';

            % Create rocket
            app.rocket = uipanel(app.UIFigure);
            app.rocket.Title = '機体諸元';
            app.rocket.Position = [10 465 479 211];

            % Create gLabel
            app.gLabel = uilabel(app.rocket);
            app.gLabel.HorizontalAlignment = 'right';
            app.gLabel.Interpreter = 'tex';
            app.gLabel.Position = [256 134 64 22];
            app.gLabel.Text = '固定質量';

            % Create constant_weight
            app.constant_weight = uieditfield(app.rocket, 'numeric');
            app.constant_weight.ValueChangedFcn = createCallbackFcn(app, @constant_weightValueChanged, true);
            app.constant_weight.Position = [335 133 60 25];
            app.constant_weight.Value = 26.3;

            % Create NsecLabel
            app.NsecLabel = uilabel(app.rocket);
            app.NsecLabel.HorizontalAlignment = 'right';
            app.NsecLabel.Interpreter = 'tex';
            app.NsecLabel.Position = [10 24 173 22];
            app.NsecLabel.Text = 'トータルインパルス/N*sec';

            % Create total_impulse
            app.total_impulse = uieditfield(app.rocket, 'numeric');
            app.total_impulse.ValueChangedFcn = createCallbackFcn(app, @total_impulseValueChanged, true);
            app.total_impulse.Position = [186 21 41 28];
            app.total_impulse.Value = 10;

            % Create sEditFieldLabel_2
            app.sEditFieldLabel_2 = uilabel(app.rocket);
            app.sEditFieldLabel_2.HorizontalAlignment = 'right';
            app.sEditFieldLabel_2.Interpreter = 'tex';
            app.sEditFieldLabel_2.Position = [91 63 62 22];
            app.sEditFieldLabel_2.Text = '噴射時間/s';

            % Create thrust_time
            app.thrust_time = uieditfield(app.rocket, 'numeric');
            app.thrust_time.ValueChangedFcn = createCallbackFcn(app, @thrust_timeValueChanged, true);
            app.thrust_time.Position = [175 61 51 27];
            app.thrust_time.Value = 1.6;

            % Create Label_2
            app.Label_2 = uilabel(app.rocket);
            app.Label_2.HorizontalAlignment = 'right';
            app.Label_2.Interpreter = 'tex';
            app.Label_2.Position = [23 101 130 34];
            app.Label_2.Text = '機体の抵抗係数(Cd)';

            % Create drag_coefficient
            app.drag_coefficient = uieditfield(app.rocket, 'numeric');
            app.drag_coefficient.ValueChangedFcn = createCallbackFcn(app, @drag_coefficientValueChanged, true);
            app.drag_coefficient.Position = [157 108 69 27];
            app.drag_coefficient.Value = 0.3;

            % Create Label_4
            app.Label_4 = uilabel(app.rocket);
            app.Label_4.HorizontalAlignment = 'right';
            app.Label_4.Interpreter = 'tex';
            app.Label_4.Position = [23 144 98 22];
            app.Label_4.Text = '機体断面積/m^{2}';

            % Create S_vertical
            app.S_vertical = uieditfield(app.rocket, 'numeric');
            app.S_vertical.ValueChangedFcn = createCallbackFcn(app, @S_verticalValueChanged, true);
            app.S_vertical.Position = [157 143 69 28];
            app.S_vertical.Value = 0.000625;

            % Create kgLabel
            app.kgLabel = uilabel(app.rocket);
            app.kgLabel.HorizontalAlignment = 'right';
            app.kgLabel.Interpreter = 'tex';
            app.kgLabel.Position = [225 24 96 22];
            app.kgLabel.Text = '燃焼後質量';

            % Create fuel_weight_end
            app.fuel_weight_end = uieditfield(app.rocket, 'numeric');
            app.fuel_weight_end.ValueChangedFcn = createCallbackFcn(app, @fuel_weight_endValueChanged2, true);
            app.fuel_weight_end.Position = [342 21 53 28];
            app.fuel_weight_end.Value = 12.48;

            % Create Label_11
            app.Label_11 = uilabel(app.rocket);
            app.Label_11.HorizontalAlignment = 'right';
            app.Label_11.Interpreter = 'tex';
            app.Label_11.Position = [226 63 96 22];
            app.Label_11.Text = '点火時質量';

            % Create fuel_weight_begin
            app.fuel_weight_begin = uieditfield(app.rocket, 'numeric');
            app.fuel_weight_begin.ValueChangedFcn = createCallbackFcn(app, @fuel_weight_beginValueChanged, true);
            app.fuel_weight_begin.Position = [342 60 53 25];
            app.fuel_weight_begin.Value = 25.8;

            % Create checking_kg_or_g
            app.checking_kg_or_g = uidropdown(app.rocket);
            app.checking_kg_or_g.Items = {'kg', 'g'};
            app.checking_kg_or_g.ValueChangedFcn = createCallbackFcn(app, @checking_kg_or_gValueChanged, true);
            app.checking_kg_or_g.Position = [415 127 54 36];
            app.checking_kg_or_g.Value = 'g';

            % Create Panel_launch_deg
            app.Panel_launch_deg = uipanel(app.UIFigure);
            app.Panel_launch_deg.Title = '打ち上げ設定';
            app.Panel_launch_deg.Position = [10 21 365 175];

            % Create ZdegLabel
            app.ZdegLabel = uilabel(app.Panel_launch_deg);
            app.ZdegLabel.HorizontalAlignment = 'right';
            app.ZdegLabel.Interpreter = 'tex';
            app.ZdegLabel.Position = [11 99 123 22];
            app.ZdegLabel.Text = '打ち上げ角度\theta/deg';

            % Create theta
            app.theta = uieditfield(app.Panel_launch_deg, 'numeric');
            app.theta.ValueChangedFcn = createCallbackFcn(app, @thetaValueChanged, true);
            app.theta.Position = [140 93 73 34];

            % Create phidegEditFieldLabel
            app.phidegEditFieldLabel = uilabel(app.Panel_launch_deg);
            app.phidegEditFieldLabel.HorizontalAlignment = 'right';
            app.phidegEditFieldLabel.Interpreter = 'tex';
            app.phidegEditFieldLabel.Position = [10 58 125 22];
            app.phidegEditFieldLabel.Text = '打ち上げ角度\phi/deg';

            % Create phi
            app.phi = uieditfield(app.Panel_launch_deg, 'numeric');
            app.phi.ValueChangedFcn = createCallbackFcn(app, @phiValueChanged, true);
            app.phi.Position = [141 52 73 34];

            % Create Image
            app.Image = uiimage(app.Panel_launch_deg);
            app.Image.ScaleMethod = 'stretch';
            app.Image.Position = [225 17 123 110];
            app.Image.ImageSource = fullfile(pathToMLAPP, '1235_0_rthetaphi1.png');

            % Create Button_2
            app.Button_2 = uibutton(app.UIFigure, 'push');
            app.Button_2.ButtonPushedFcn = createCallbackFcn(app, @Button_2Pushed, true);
            app.Button_2.Position = [945 225 129 49];
            app.Button_2.Text = 'グラフの初期化';

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