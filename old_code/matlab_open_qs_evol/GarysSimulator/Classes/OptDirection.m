classdef OptDirection < handle
    properties
        oldDirection
        newDirection
        
        oldGradient
        newGradient
        
        deltaPulse %correspond to variation in pulse of previous update
        
        oldHessianInv
        newHessianInv
        
    end
    
    methods
        function obj = OptDirection()
            obj.oldDirection = 0;
            obj.newDirection = 0;

            obj.oldGradient = 0;
            obj.newGradient = 0;

            obj.oldHessianInv = 0;
            obj.newHessianInv = 0;
        end
        
        function directionImprov = Next(obj, type)
            %Define DIRECTION of improvement
            
            HesLength = numel(obj.newGradient);
                    
            %Direction for first iteration
            if(obj.oldGradient == 0)
                obj.newDirection = obj.newGradient;
                if(strcmp(type,'BFGS')) 
                    obj.oldHessianInv = -eye(HesLength);
                end
            else                
                switch(type)
                    case 'Gradient'
                        obj.newDirection = obj.newGradient;

                    case 'Conjugate'
                        %Polak-Ribiere method, auto direction reset
                        Dgrad = obj.newGradient - obj.oldGradient;
                        
                        beta = (obj.newGradient(:).'*Dgrad(:))/(obj.oldGradient(:).'*obj.oldGradient(:));
                        beta = max(beta,0);
                        %beta must be positive, if negative than direction becomes gradient

                        obj.newDirection = beta*obj.oldDirection + obj.newGradient;
                        
                    case 'BFGS'
                        Dgrad = obj.newGradient - obj.oldGradient;
                        divGrad = obj.deltaPulse(:).'*Dgrad(:);
                    
                        if(divGrad < 0) %Definite negative condition
                            QHes = eye(HesLength) - (Dgrad(:)*obj.deltaPulse(:).')/divGrad;
                            obj.newHessianInv = QHes.'*obj.oldHessianInv*QHes + (obj.deltaPulse(:)*obj.deltaPulse(:).')/divGrad;

                            obj.newDirection(:) = -obj.newHessianInv*obj.newGradient(:);

                            obj.oldHessianInv = obj.newHessianInv;

                        else %If divide by zero use simple gradient
                            obj.newDirection = obj.newGradient;
                            obj.oldHessianInv = -eye(HesLength); 
                        end
                        
                        %direction tends to increase close to convergence.
                        %When equal inf, convergence is achieved
                        if(isinf(max(obj.newDirection(:))))
                            obj.oldDirection = 0;
                            obj.newDirection = 0;
                        end
                        
                    otherwise
                        ERROR('METHOD NAME IS UNKNOWN');
                end
            end
                       
            %Improv direction is for algorithm stop criteria
            meanOld = mean(abs(obj.oldDirection(:)));
            meanNew = mean(abs(obj.newDirection(:)));
            if(meanOld == meanNew)
                directionImprov = 0;
            elseif(meanOld ~= 0)
                directionImprov = abs(1 - meanNew/meanOld);
            else
                directionImprov = 1;
            end
            
            %Memorize old direction
            obj.oldGradient = obj.newGradient;
            obj.oldDirection = obj.newDirection;
                 
        end
        
    end
    
end

