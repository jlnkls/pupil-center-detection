function [ H ] = ILPF(size, D0)
%ILPF Two-dimensional ideal low-pass filter (LPF)
% 
%   Input:
%       size                        filter size
%       D0                          circumference radius of the passband
%
%   Output:
%       H                           filter (final matrix)
%
% 
% 
%   Author: jlnkls
%
% 	26/01/2016


%% Transfer function
[x,y] = meshgrid(1:size(1), 1:size(2));   

H = sqrt( ((x-(size(1)/2)).^2) + ((y-(size(2)/2)).^2) );

%% Filter implementation
H((H<D0)) = 1;
H((H>=D0)) = 0;

end
