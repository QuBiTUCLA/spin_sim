classdef PlotSpins < handle
    %Allow "measurements" of the system by propagation and projection at
    %each time step for the spin watched
    
    properties
        Operators
        Results
        TotalPulse
        Normalize
        WatchNames
    end
    
    methods
        
        function obj = PlotSpins()
        end
        
        %Will measure and plot the operator
        function WatchSpin(obj, Name, Hnat)
            endIndex = length(obj.Operators)+1;
            spin = Hnat.findObject(Name);
            spinDim = spin.dimH;
            obj.WatchNames{endIndex} = Name;
            
            %No normalization as done in Propagate method
            switch(spinDim)
                case 2
                    %ms=0
                    obj.Operators{endIndex}{1} = Hnat.expandOperator(Name,diag([0 1]));

                    %ms=1                                                  
                    obj.Operators{endIndex}{2} = Hnat.expandOperator(Name,diag([1 0]));

                case 3
                    %ms=0
                    obj.Operators{endIndex}{1} = Hnat.expandOperator(Name,diag([0 1 0]));

                    %ms=1                                                  
                    obj.Operators{endIndex}{2} = Hnat.expandOperator(Name,diag([1 0 0]));

                    %ms=-1                                                  
                    obj.Operators{endIndex}{3} = Hnat.expandOperator(Name,diag([0 0 1]));    
                otherwise
                    ERROR('Can only watch spin 2 and 3');
            end
            
            %ExpandOperator is equal to get a mixed state for the environment
            %So trace will be equivalent to partial trace                                                           
        end
        
        function Propagate(obj, propagator)
            %propagator is the simulator class with all control matrices
            %pulse is the control pulse
            
            rhoin = propagator.rhoIn;
            timeLength = propagator.TimeLength;
            
            %MANY GOAL NOT IMPLEMENTED
            obj.TotalPulse = (0:timeLength)*propagator.TimeStep;
            
            obj.Results = cell(size(obj.Operators));
            obj.Normalize = zeros(size(obj.Operators));
            
            for endIndex = 1:length(obj.Operators);
                obj.Results{endIndex} = cell(size(obj.Operators{endIndex}));
            %Results: projection on different ms states
            %1st step
                for watchIndex = 1:length(obj.Operators{endIndex})
                    obj.Results{endIndex}{watchIndex}(1) = abs(trace(obj.Operators{endIndex}{watchIndex}'*rhoin));
                    obj.Normalize(endIndex) = obj.Normalize(endIndex) + obj.Results{endIndex}{watchIndex}(1);
                end
            end
            
            %Prepare Final Hamiltonian
            U = propagator.unitary();
            rhoout = rhoin;

            %Propagate
            j = 2; %rho0 -> rhoN = N+1 in the end
            for ct=1:timeLength
                rhoout = U{ct}*rhoout*U{ct}';
                for endIndex = 1:length(obj.Operators)
                %Results: projection on different ms states
                    for watchIndex = 1:length(obj.Operators{endIndex})
                        obj.Results{endIndex}{watchIndex}(j) = abs(trace(obj.Operators{endIndex}{watchIndex}'*rhoout));
                    end
                end
                j = j + 1;
            end       
        end
        
        function Plot(obj)
%             clf
            for endIndex=1:length(obj.Operators)
                figure(endIndex);
                
                hold on
                plot(obj.TotalPulse*1e6,obj.Results{endIndex}{1}/obj.Normalize(endIndex),'*-r');
                plot(obj.TotalPulse*1e6,obj.Results{endIndex}{2}/obj.Normalize(endIndex),'*-b');
                if(length(obj.Results{endIndex}) == 3)
                    plot(obj.TotalPulse*1e6,obj.Results{endIndex}{3}/obj.Normalize(endIndex),'*-g');
                end
                axis([0 obj.TotalPulse(1,end)*1e6 0 1]);
                hold off

                xlabel('Pulse duration (us)');
                text = strcat('Probabilities for',' ',obj.WatchNames{endIndex},' states');
                ylabel(text);
                legend('m=0','m=1')
            end
            
        end
        
        function FFTPlot(obj) %Change if the duration is now in s and not in rad
                figure(3);
                  
                L = size(obj.Results{1}{1},2);
                Fs = 1/(obj.TotalPulse(2)-obj.TotalPulse(1));%Sampling frequency
                
                NFFT = 2^nextpow2(L); % Next power of 2 from length of y
                Y = fft(obj.Results{1}{1}-mean(obj.Results{1}{1}),NFFT)/L;
                f = Fs/2*linspace(0,1,NFFT/2+1);

                % Plot single-sided amplitude spectrum.
                plot(f,2*abs(Y(1:NFFT/2+1))) 

                xlabel('Pulse frequencies(hz)');
                ylabel('FFT of m=0 Rabi Oscillation');
                legend('m=0')
        end
                
    end
    
end