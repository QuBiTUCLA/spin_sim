classdef PulseSimObj
    %Can be used for free precession by defining RotFrame
    
    properties
        Pulse
        Size %Control Field and Time size
        RotFrame
    end
    
    methods
        
        function obj = PulseSimObj(sizeFieldTime)
            obj.Size = sizeFieldTime;
            obj.Pulse = zeros(sizeFieldTime);
            obj.RotFrame = ones(sizeFieldTime(2),1); %1 normal, 2 for RWA
        end
        
         %Create a random pulse 
        function outPulse = makeRandomPulse(obj,randPower,pulseScale)
            %Create a skeleton of random points and then
            %fit the points with a interpolating cubic spline

            if(randPower > 1)
               error('randPower must be <1');
            end
            randLength = ceil(obj.Size(2)*randPower)+1; %2 data point minimum
            randDuration = linspace(1,obj.Size(2),randLength);
            pulseDuration = 1:obj.Size(2);
            
            outPulse = zeros(obj.Size(1),obj.Size(2));
            for ctField=1:obj.Size(1)
                randPulse = pulseScale*(2*rand(randLength,1)-1);
                outPulse(ctField,:) = interp1(randDuration,randPulse,pulseDuration,'spline');
            end
        end
    end   
end