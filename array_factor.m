function [AF,theta,phi] = array_factor(r_prime,phase_correction,resolution)
%array_factor Summary of this function goes here
%INPUT x(n) = sin(2piny*n + arbitrary phase)     
%OUTPUT y(n) = sin(2piny*n + arbitrary phase + phase)
%   Detailed explanation goes here
theta = linspace(0,pi/2,resolution);
phi = linspace(0,2*pi,resolution);
AF = zeros(length(theta),length(phi));

for i = 1:length(theta)
    for j = 1:length(phi)
        sum = 0;
        r_hat = r_vec(theta(i),phi(j));
        for element = 1:length(r_prime(1,:))
            sum = sum + exp(1i*(dot(r_hat,r_prime(:,element)))) * ...
                phase_correction(element);
        end
        AF(i,j) = sum;
    end
end