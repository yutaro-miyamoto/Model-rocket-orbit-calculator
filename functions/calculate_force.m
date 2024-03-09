function output = calculate_force(weight,gravitational_acceleration,air_resistance, Propellant_Force)
output = [0;0;-weight *gravitational_acceleration] -air_resistance + Propellant_Force;
end

