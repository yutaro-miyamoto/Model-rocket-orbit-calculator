function output = velocity(verocity_before, average_exhaust_velocity, Weight_of_fuel_using, Weight, force, dt)
output =  verocity_before + average_exhaust_velocity * Weight_of_fuel_using / Weight +  force / Weight * dt;
end