function output_abs = GSS_MPDR_angle_obj_func(ang, frequency, freqs_vector, c, MicPos, Ryy, MicNum, dia_load_beamformer)

output_abs = 0;
kappa = [sind(ang), cosd(ang), 0];
for n = 1:frequency
    omega = 2*pi*freqs_vector(n);
    steer_vec = exp(1j*omega/c*kappa*MicPos.').';
    array_output_power = 1/(steer_vec'*inv(Ryy(:, :, n)+dia_load_beamformer*eye(MicNum))*steer_vec);
    output_abs = output_abs + abs(array_output_power);

end
