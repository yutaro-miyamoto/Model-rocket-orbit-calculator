classdef calculation
    %AGEO このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties%データを用意
        velocity
        location
        mag
        heading_vector
        air_resistance
    end

    properties (Access = private)
        force
    end

    methods(Static)

        function obj = ...
                calculation(Weight, gravitational_acceleration, air_resistance, velocity_before, ...
                location_before, Weight_of_fuel_using, wind_speed, dt, Isp, heading_vector, drag_coefficient,rho,cross_sectional_area)
            %WHEN_ROCKET_ON この関数の概要をここに記述
            %ロケット推進がOnの時の弾道計算を行うメソッド
            force = [0;0;-Weight * gravitational_acceleration] - air_resistance;
            % 速度
            obj.velocity = velocity_before +  force / Weight * dt + Isp * gravitational_acceleration * log(Weight / (Weight-Weight_of_fuel_using)) * heading_vector+ wind_speed;
            % 位置
            obj.location = location_before + obj.velocity * dt;
            % 速度絶対値
            obj.mag = sqrt(sum(obj.velocity.^2));
            % 機首方向単位ベクトル
            obj.heading_vector =  obj.velocity / obj.mag;
            % 空気抵抗
            obj.air_resistance = 0.5 * drag_coefficient * rho * sqrt(sum(obj.velocity.^2)) * obj.velocity * cross_sectional_area;
        end
    end
end

