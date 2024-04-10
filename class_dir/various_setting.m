classdef various_setting
    %VARIOUS_SETTING このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        m_fuel
        m
        m_fuel_using
        F_r0
        Isp
    end

    methods
        function obj =...
                various_setting(fuel_weight_begin, fuel_weight_end, times_const,constant_weight,...
                thrust_time, total_impulse, dt)
            %　総消費燃料質量(消費される推進剤の全質量)
            obj.m_fuel = fuel_weight_begin - fuel_weight_end;
            obj.m_fuel = obj.m_fuel * times_const;%kgへの単位変換

            % 全質量/g
            obj.m = constant_weight + fuel_weight_begin;
            obj.m = obj.m * times_const;%kgへの単位変換

            %(dt)s当たりの消費燃料質量
            obj.m_fuel_using = obj.m_fuel / thrust_time * dt;

            %推力/N
            obj.F_r0 = [0;0;total_impulse / thrust_time];

            %平均比推力Isp
            obj.Isp = total_impulse / obj.m_fuel;
        end

    end
end

