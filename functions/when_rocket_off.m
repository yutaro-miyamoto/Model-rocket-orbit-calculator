function [velocity,location] = ...
    when_rocket_off(Weight, gravitational_acceleration, air_resistance, velocity_before, ...
                location_before,wind_speed, dt)
%WHEN_ROCKET_OFF この関数の概要をここに記述
%ロケット推進がoffのときの軌道計算はここをいじる。
%   詳細説明をここに記述
force = [0;0;-Weight *gravitational_acceleration] - air_resistance;

velocity = velocity_before +   force / Weight * dt;

location = location_before + (velocity + wind_speed) * dt;
end

