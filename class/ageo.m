classdef ageo
    %AGEO このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties%データを用意
        velocity;
        location;
    end

    properties (Access = private)
        force;
    end

    methods(Static)
        function obj = ...
                when_rocket_on(Weight, gravitational_acceleration, air_resistance, velocity_before, ...
                location_before, Weight_of_fuel_using, wind_speed, dt, Isp, heading_vector)
            %WHEN_ROCKET_ON この関数の概要をここに記述
            %ロケット推進がOnの時の弾道計算を行う。
            %   詳細説明をここに記述
            obj.force = [0;0;-Weight * gravitational_acceleration] - air_resistance;

            obj.velocity = velocity_before +  obj.force / Weight * dt + Isp * gravitational_acceleration * log(Weight/(Weight-Weight_of_fuel_using)) * heading_vector + wind_speed;

            obj.location = location_before + obj.velocity * dt;
        end
    end
end

