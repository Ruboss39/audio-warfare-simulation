classdef audio_source
    % Audio Source
    %   Detailed explanation goes here
    
    properties
        theta
        phi
        rho
        frequency
        t_start
        t_end

    end
    
    methods
        function obj = audio_source(f_start,f_end,f_res,theta_deg,phi_deg,rho,t_start,t_end)
            %MATRIX_ARRAY Construct an instance of this class
            %   Detailed explanation goes here
            %
            %   Audio source with directional position data
            %   Audio source is a sum of sine waves with frequencies
            %   ranging from f_start to f_end, with uniform distribution
            %   
            %   Rho is the distance between the origin and the audio source
            %
            %   t_start and t_end is the span of time of which the audio
            %   source emits sound
            %
            %
            obj.theta = theta_deg*pi/180;
            obj.phi = phi_deg*pi/180;
            obj.frequency = linspace(f_start,f_end,f_res);
            obj.t_start = t_start;
            obj.t_end = t_end;
            obj.rho = rho;
        end
    end
end

