classdef various_setting
    %VARIOUS_SETTING このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        weight_fuel
        weight
        weight_fuel_using
        thrust_power
        Isp
    end

    methods
        function obj =...
                various_setting(fuel_weight_begin, fuel_weight_end, times_const,constant_weight,...
                thrust_time, total_impulse, dt, gravitational_acceleration)
            
            %　総消費燃料質量(消費される推進剤の全質量)
            obj.weight_fuel = fuel_weight_begin - fuel_weight_end;
            obj.weight_fuel = obj.weight_fuel * times_const;%kgへの単位変換

            % 全質量/g
            obj.weight = constant_weight + fuel_weight_begin;
            obj.weight = obj.weight * times_const;%kgへの単位変換

            %(dt)s当たりの消費燃料質量
            obj.weight_fuel_using = obj.weight_fuel / thrust_time * dt;

            %推力/N
            obj.thrust_power = [0;0;total_impulse / thrust_time];
            
            %平均比推力Isp/s
            obj.Isp = (total_impulse/ thrust_time) / (obj.weight_fuel * gravitational_acceleration) ;%(kgで処理しなくてはいけない)
        end

    end
end

