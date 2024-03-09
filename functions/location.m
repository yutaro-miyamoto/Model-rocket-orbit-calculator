function output = location(location_before, velocity, wind_speed, dt)
output = location_before + (velocity + wind_speed) * dt;
end

