function [ y ] = simpson_solver( f, h, L,L_cap, V, height, lat, lon, t, mu_r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    x = 0:h:L;
    y = 0;
    
    for i = 1:length(x)-1
        xi = x(i);
        y = y + h*(f(L,L_cap, V, height - xi, lat, lon, t, mu_r) + ...
            4*f(L,L_cap, V, height - xi + h/2, lat, lon, t, mu_r) + ...
            f(L,L_cap, V, height - xi + h, lat, lon, t, mu_r))/6;
    end
    
end

