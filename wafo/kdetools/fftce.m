function y = fftce(x)
% FFTCE  Circulant Embedding of a vector or matrix
%
% CALL: y = fftce(x);
%
%   y = circulant embedded vector or matrix (2*size(x)-1)
%   x = vector or matrix 
%
%   For vectors, Y = X([1:M M-1:-1:2])
%   For matrices, X is copied into the first quadrant of Y
%   and Mi-1 on to the 2'nd, 3'rd and 4'th quadrant of Y.
%   For N-D arrays, X is copied into  "half-spaces" of Y along each
%   dimension.
%
%   FFTCE is useful for circulant embedding of vectors and matrices
%   when working with the Fourier transform. 
%
%   See also  fftshift, ifftshift, fft, fft2, fftn

% tested on: matlab 5.2
% history:
% by pab 5.11.1999

numDims = ndims(x);
idx = cell(1, numDims);
for k = 1:numDims
    m = size(x, k);
    idx{k} = [1:m m-1:-1:2];
    %idx{k} = [1:m m:-1:2];
end

% Use comma-separated list syntax for N-D indexing.
y = x(idx{:});
