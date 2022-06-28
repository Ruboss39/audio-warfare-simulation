function [x,y,z,x_dbi,y_dbi,z_dbi] = farfield(E_mag,theta,phi,sensitivity)
%array_factor Summary of this function goes here
%INPUT: magnitude of the E-field (matrix theta x phi), theta vector and phi
%vector. Theta and phi vector are column vectors. Sensitivity given in dB, 
%OUTPUT: x y z coordinates for 3d plotting
%   Detailed explanation goes here
x = zeros(length(theta),length(phi));
y = zeros(length(theta),length(phi));
z = zeros(length(theta),length(phi));

power = abs(E_mag).^2;
integrand = ((sin(theta))').*power;

temp = cumtrapz(phi,cumtrapz(theta,integrand',2));
power_rad = (temp(end,end) - temp(1,1))/(4*pi);

E_mag_dbi = 10*log10(power/power_rad); %Normalized power

for i = 1:length(theta)
    for j = 1:length(phi)
        x(i,j) = abs(E_mag(i,j)) * (sin(theta(i))*cos(phi(j)));
        y(i,j) = abs(E_mag(i,j)) * (sin(theta(i))*sin(phi(j)));
        z(i,j) = abs(E_mag(i,j)) * cos(theta(i));
    end
end

for i = 1:length(E_mag_dbi(:,1))
    for j = 1 :length(E_mag_dbi(1,:))

        E_mag_dbi(i,j) = E_mag_dbi(i,j) - sensitivity;
        if(E_mag_dbi(i,j) < 0)
            E_mag_dbi(i,j) = 0;
        end
    end
end

x_dbi = zeros(length(theta),length(phi));
y_dbi = zeros(length(theta),length(phi));
z_dbi = zeros(length(theta),length(phi));

for i = 1:length(theta)
    for j = 1:length(phi)
        x_dbi(i,j) = abs(E_mag_dbi(i,j)) * (sin(theta(i))*cos(phi(j)));
        y_dbi(i,j) = abs(E_mag_dbi(i,j)) * (sin(theta(i))*sin(phi(j)));
        z_dbi(i,j) = abs(E_mag_dbi(i,j)) * cos(theta(i));
    end
end


