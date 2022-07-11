function mic_data = beam_forming_alogrithm(matrix_array,direction,weight,audio_signal,frequency,sampling_frequency,c)
%   
%   IMPORTANT!!!! The beamforming algorithm assumes the array matrices lies in the
%   xy-plane
%
%   The beamforming algorithm calculates the necessary phase to introduce
%   to the narrowband signal in order to have a maximum directivity in the
%   direction r(theta,phi)
%
%   This phase shift is introduced by a phase shifting function, which acts
%   as a filter in time domain
%
%   To improve performance, all elements of the array matrices are not in
%   use. The user decides which element to use by sending in a weight
%   vector as an argument. The output signal is also then normalized after
%   how many elements where in use.
%

% Get the amount of samples of the audio tracks
samples = length(audio_signal(:,1));

% The listening direction vector contains two scalar values, theta and phi
theta = direction(1);
phi = direction(2);

% The r_prime vector of the matrix array to know the location of every
% element, as well as how many elements exists.
r_prime = matrix_array.r_prime;
elements = matrix_array.elements;

% The narrowband frequency
k = 2*pi*frequency/c;

% The normalized frequency
ny = frequency/sampling_frequency;

% Initiazte output vector
mic_data = zeros(samples,1);

% The compensation factors to obtain uniform phase in the direction
% r_hat(theta,phi)
x_factor = sin(theta)*cos(phi);
y_factor = sin(theta)*sin(phi);

for mic_ind = 1:elements

    % Calculate the narrowband phase-shift
    phase_shift_value = -k*(r_prime(1,mic_ind)*x_factor + ...
         r_prime(2,mic_ind)*y_factor);

    %phase_shift_value = phase_shift_value - floor(phase_shift_value/(2*pi));

    % Sum the individually shifted data from the antenna elements, as well
    % as weight them with the appropriate weight.
    mic_data = mic_data + weight(mic_ind)*phase_shift(audio_signal(:,mic_ind),ny,phase_shift_value);
end

% Normalized output
norm_coeff = 1/sum(weight);
mic_data = mic_data*norm_coeff;