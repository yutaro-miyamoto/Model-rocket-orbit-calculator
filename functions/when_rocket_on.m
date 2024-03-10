function [velocity, location] = ...
    when_rocket_on(Weight, gravitational_acceleration, air_resistance, Propellant_Force, velocity_before, ...
                location_before, Weight_of_fuel_using, average_exhaust_velocity,wind_speed, dt)
%WHEN_ROCKET_ON この関数の概要をここに記述
%ロケット推進がOnの時の弾道計算を行う。
%   詳細説明をここに記述
force = [0;0;-Weight *gravitational_acceleration] - air_resistance + Propellant_Force;

velocity = velocity_before + average_exhaust_velocity * Weight_of_fuel_using / Weight +  force / Weight * dt;

location = location_before + (velocity + wind_speed) * dt;
end

