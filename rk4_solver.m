function [ y ] = rk4_solver( f,x0, xf, h, y0 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    x = x0:h:xf;
    y = y0;
    
    for i = 1:length(x)-1
        xi = x(i);
        k1 = h*f(xi,y);
        k2 = h*f(xi + h/2, y+k1/2); 
        k3 = h*f(xi + h/2, y+k2/2);
        k4 = h*f(xi + h, y+k3);
        y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

end

