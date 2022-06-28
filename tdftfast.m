function [Y, nu] = tdftfast(y)
Y = fft(y);
N = length(y);
nu = [(0:N-1)./N];
Y = fftshift(Y);
nu = fftshift(nu);
nu_index = find(nu >= 0.5);
nu(nu_index) = nu(nu_index)-1;