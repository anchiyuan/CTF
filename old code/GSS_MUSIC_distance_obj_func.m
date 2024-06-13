function output_power_sum = GSS_MUSIC_distance_obj_func(dis, ang, MicNum, MicPos, frequency_lower_bound, frequency, freqs_vector, c, P_N)

output_power_sum = 0;
source_pos = [dis*sind(ang), dis*cosd(ang), 0];
r = zeros(MicNum, 1);
for i = 1 : MicNum
    r(i, :) =  sqrt(sum((source_pos - MicPos(i, :)).^2));

end

for n = frequency_lower_bound:frequency
    omega = 2*pi*freqs_vector(n);
    steer_vec = exp(-1j*omega/c*r)./r;

    array_output_power = 1/(steer_vec'*P_N(:, :, n)*steer_vec);
    output_power_sum = output_power_sum + abs(array_output_power);
    
end
