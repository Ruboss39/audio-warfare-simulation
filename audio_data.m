classdef audio_data
    %AUDIO_DATA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        audio_signals
    end
    
    methods
        function obj = audio_data(audio_signals)
            %AUDIO_DATA Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0
            obj.audio_signals = audio_signals;
            end
            
        end
    end
end

