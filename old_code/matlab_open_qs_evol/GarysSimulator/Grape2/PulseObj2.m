classdef PulseObj2
    
    properties
        Pulse
        Index
        Size
    end
    
    methods
        
        function obj = PulseObj2(sizeFieldTime)
            obj.Pulse = zeros(sizeFieldTime);
            obj.Index = zeros(sizeFieldTime);
            obj.Size = sizeFieldTime;
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
        
        %Create a random pulse index (INTEGER VALUES)
        function outPulse = makeRandomIndex(obj,randPower,pulseScale)
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
                randPulse = pulseScale(1)+(pulseScale(2)-pulseScale(1))*rand(randLength,1);
                outPulse(ctField,:) = round(interp1(randDuration,randPulse,pulseDuration,'spline'));
                outPulse(ctField,:) = max(pulseScale(1),min(pulseScale(2),outPulse(ctField,:)));
            end
        end
    end   
end