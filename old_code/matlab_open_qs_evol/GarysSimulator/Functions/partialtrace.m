function [output] = partialtrace(input,dim)
%Give the partial trace of a matrix. Output is thus a matrix
%Input: matrix to trace and dimension of partial trace

%Must be normalized!!!!
%VERY LIMITED

if(size(input,1) ~= size(input,2))
    ERROR('Matrix input must be square')
end
dimEnv = size(input,1)/dim;
if(dimEnv ~= round(dimEnv))
    ERROR('Dimension of input matrix must be a multiple of dim')
end

output = zeros(dim);

for k=1:dim
    for i=1:dim
        output(k,i) = trace(input(1+(k-1)*dimEnv:k*dimEnv,1+(i-1)*dimEnv:i*dimEnv));
    end
end

end

