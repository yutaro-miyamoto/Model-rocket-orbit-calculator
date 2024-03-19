function air_resistance = air_resistance(drag_coefficient, rho, velocity, cross_sectional_area)
% 空気抵抗
air_resistance = 0.5 * drag_coefficient * rho * sqrt(sum(velocity.^2)) * velocity * cross_sectional_area;

end

