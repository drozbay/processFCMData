function [ S ] = roundOdd( S )
%ROUNDODD Summary of this function goes here
%   Detailed explanation goes here
idx = mod(S,2)<1;
S = floor(S);
S(idx) = S(idx)+1;

end

