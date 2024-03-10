function [heading_vector, mag] = heading_vector_fuction(velocity)
% 速度ベクトルを用いたヘディングベクトルの生成
mag = sqrt(sum(velocity.^2));
heading_vector =  velocity / mag;
end

