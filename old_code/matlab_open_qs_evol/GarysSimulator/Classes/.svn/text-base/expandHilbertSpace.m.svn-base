function [H2] = expandHilbertSpace(H,newLogicalIndicies,newSpaceIndicies,dimensions)
%function [H2] = expandHilbertSpace(H,newLogicalIndicies,newSpaceIndicies,dimensions)
%
% 31May2006
% jhodges
% Take a smaller matrix and creates a larger matrix, assuming identity on
% other spins in the larger hilbert space
% 
% 12Oct2008
% Revised to deal with subsystems of arbitrary size
%
% 5Mar2010
% Colm Ryan - minor cleaning of unused variables
%
% Usage:
% H: smaller matrix
% 
%newLogicalIndicies:   give the ordering in the H matrix, how to reorder
% in the H2 matrix
%
%newSpaceIndicies:  the parts of the new space which are in identity
%
%dimensions: (optional) specify the dimension for each subsystem as a
%single vector
%
%
% Example:
%
% Suppose you have a 4x4 matrix that represents the Hilbert space between
% spins 1&3.  You want to represent the matrix on the 8x8 Hilbert space
% between spins 1,2, & 3 but you cannot use the kron() function because the
% ordering of the spins is not reflected in the 4x4 matrix
%
% use H8 = expandHilbertSpace(H4,[1,3],[2]) to get the right answer.
%
% Note expandHilbertSpace(H4,[1,2],[3]) = kron(H4,eye(2));
%
% Technical Note:
%
% kron(A,B) is not equal to the kron(B,A)
% However, the two are related by a (symmetric) permutation matrix
% In expanding the Hilbert space, we calculate this permutation matrix and
% apply to a simple kron(H,I)


% if no dimensions given, assume spin-1/2 or qubit
if nargin < 4,
    dimensions = 2*ones(length([newLogicalIndicies,newSpaceIndicies]),1);
end


% matrix size
hSize = prod(dimensions(newLogicalIndicies));

% new matrix size
h2Size = prod(dimensions);

% make the new matrix in the wrong basis
H2 = kron(H,eye(h2Size/hSize));

% make ordered list of indicies
OldIndicies = [newLogicalIndicies,newSpaceIndicies];

% T is the total transformation matrix which reorders the kron products
T = eye(h2Size);

% loop over the Indicies
for k=1:length(OldIndicies),
    % if not ordered correctly, swap the order pairwise until the end
    while  k < OldIndicies(k),
        
        % find where the kth entry is
        curInd = find(OldIndicies == k);
        
        % calculate the permutation matrix for swaping the matricies of the
        % from kron(I,B,A,I) to kron(I,A,B,I)
        dimA = dimensions(OldIndicies(curInd));
        dimB = dimensions(OldIndicies(curInd-1));
        P = CalcPermMatrix([dimA,dimB]);
        PreDim = prod(dimensions(OldIndicies(1:curInd-2)));
        PostDim = prod(dimensions(OldIndicies(curInd+1:end)));
        Pre = eye(max(1,PreDim));
        Post = eye(max(1,PostDim));
        T = kron(Pre,kron(P,Post))*T;
        
        % reorder OldIndicies to reflect the change
        OldIndicies([curInd-1,curInd]) = OldIndicies([curInd,curInd-1]);
    end
end

%Apply the permutation matrix
H2 = T*H2*T';

function [P] = CalcPermMatrix(dims)

P = zeros(size(sum(dims)));
for k=1:dims(1)
    for l=1:dims(2),
      P((k-1)*dims(2) + l,(l-1)*dims(1)+k) = 1;
    end
end
