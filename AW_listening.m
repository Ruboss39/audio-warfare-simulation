function audio_out = AW_listening(c,f_sampling,array_matrices,theta_listen,phi_listen,array_audio_signal)
%
%   Spatial filtering in r_hat(theta,phi)
%
%
%   This function takes in the audio_signals generated and the
%   corresponding array matrix-geometries. This function further needs some
%   information about the audio_signals such as the sampling frequency,
%   soundspeed in the medium (c).
%
%   A direction is needed (theta,phi) to output a spatially filtered
%   signal.
%
%   The function splits the input audio signals recorded by the
%   array_matrices and then perform narrowband beamforming on the signals.
%   Different audio signals corresponds to different sub_arrays.
%

% Initialize our audio out vector
samples = length(array_audio_signal(1).audio_signals(:,1));
audio_out = zeros(samples,1);

% Array data
sub_arrays = length(array_matrices);

% Convert input theta and phi in degrees to radians
theta = pi/180*theta_listen;
phi = pi/180*phi_listen;

% Which bands to check
f_bands_N = 45;          %Linearly spaced N bands,
frequency_bands = linspace(100,3900,f_bands_N);

% Create N-audio complete audio signals for every band
audio_filtered_complete(sub_arrays,f_bands_N) = audio_data;
filter_order = 200;
filter_coefficients = zeros(f_bands_N,filter_order+1);

for array = 1:sub_arrays
    % Filter the signals generated on the arrays
    audio_signal = array_audio_signal(array).audio_signals;
    elements = array_matrices(array).elements;

    for freq_ind = 1:length(frequency_bands)
        % Filter design for each band
        nu_0 = 2*frequency_bands(freq_ind)/f_sampling;
        scale_factor = 10000;
        cut_off = [nu_0 - nu_0/scale_factor, nu_0  + nu_0/scale_factor];
        b = fir1(filter_order,cut_off);
        filter_coefficients(freq_ind,:) = b;
    
        audio_temp = zeros(samples,elements); 
        for mic_ind = 1:elements
            % Apply filter on every signal recorded from the elements
            audio_temp(:,mic_ind) = filter(b,1,audio_signal(:,mic_ind));
        end
    
        % Add the complete set of element-signals
        audio_filtered_complete(array,freq_ind) = audio_data(audio_temp);
    end

end

% Apply beamforming algo. in every freq. band
for freq_ind = 1:length(frequency_bands)      
    frequency = frequency_bands(freq_ind);
    
    mic_data = 0;
    for array = 1:sub_arrays
      % Use the filtered audio signals 
      audio_temp_signals = audio_filtered_complete(array,freq_ind).audio_signals;
    
      % Adaptive configuration of the antenna array
      % Select only necessary antenna-elements to maintain small
      % beamwidth.
      % weight = adaptive_array_config_gen2(array_matrix,frequency,theta,phi,c);
      weight = adaptive_array_config(array_matrices(array),frequency,c);
      % weight = ones(1,elements);
    
      % Perform the beamforming algorithm (phase-shift input signal according to listening direction)
      mic_data = mic_data + beam_forming_alogrithm(array_matrices(array),[theta,phi],...
          weight,audio_temp_signals,frequency,f_sampling,c);
    
    end
    audio_out = audio_out + mic_data;
end

audio_out = audio_out/max(audio_out);
end
