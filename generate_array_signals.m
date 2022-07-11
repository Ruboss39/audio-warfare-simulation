function samples = generate_array_signals(matrix_array,sources,t,c)
%generate_array_signals Summary of this function goes here
%   Detailed explanation goes here
%% Generate signals on each element of the audio array
r_prime = matrix_array.r_prime;
Audio_signal = zeros(length(t),length(r_prime(1,:)));

for sample = 1:length(t)
    for mic = 1:length(r_prime(1,:))
        x_i = r_prime(1,mic);
        y_i = r_prime(2,mic);
        temp_signal_sample = 0;
        for source = 1:length(sources)
            if sources(source).t_start <= t(sample) && t(sample) < sources(source).t_end
                frequencies_ps = sources(source).frequency;
                theta_source = sources(source).theta;
                phi_source = sources(source).phi;
                rho_source = sources(source).rho;
                for freq_ind = 1:length(frequencies_ps)
                    k = 2*pi*frequencies_ps(freq_ind)/c;
                    r_1 = [x_i;y_i;0];
                    r_2 = rho_source*r_vec(theta_source,phi_source);
                    norm_factor = norm(r_2-r_1);
                    phase_offset = -k*norm_factor;
                    %phase_offset = k*(x_i*sin(theta_source)*cos(phi_source) + ...
                        %y_i*sin(theta_source)*sin(phi_source));
                    element_amplitude = 1/norm_factor;
                    %element_amplitude = 1;
        
                    temp_signal_sample = temp_signal_sample + ...
                       element_amplitude * sin(2*pi*frequencies_ps(freq_ind)*t(sample) + phase_offset);
        
                end
            end

        end
        Audio_signal(sample,mic) = temp_signal_sample;
    end
end

samples =  Audio_signal;

end
