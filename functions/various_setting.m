
function [m_fuel, m, m_fuel_using, F_r0, Isp] =...
    various_setting(fuel_weight_begin, fuel_weight_end, times_const,constant_weight,...
                    thrust_time, total_impulse, dt)
 %　総消費燃料質量(消費される推進剤の全質量)
                    m_fuel = fuel_weight_begin - fuel_weight_end;
                    m_fuel = m_fuel * times_const;%kgへの単位変換

                    % 全質量/g
                    m = constant_weight + fuel_weight_begin;
                    m = m * times_const;%kgへの単位変換

                    %(dt)s当たりの消費燃料質量
                    m_fuel_using = m_fuel / thrust_time * dt;

                    %推力/N
                    F_r0 = [0;0;total_impulse / thrust_time];

                    %平均比推力Isp
                    Isp = total_impulse / m_fuel;
end

