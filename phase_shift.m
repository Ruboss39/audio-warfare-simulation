function y = phase_shift(x,ny,phase)
%PHASE Summary of this function goes here
%INPUT x(n) = sin(2piny*n + arbitrary phase)     
%OUTPUT y(n) = sin(2piny*n + arbitrary phase + phase)
%   Detailed explanation goes here
y = zeros(length(x),1);

for i = 2:length(x)-1
    y(i) =  cos(phase)*x(i) + sin(phase)/(2*pi*ny)*(x(i+1)/2-x(i-1)/2);
end


