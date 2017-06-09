function [ y ] = ecif2ecef( x,t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    v = [x(1) x(2)]';
    global w_earth;
    theta = w_earth*t;
    R = [cos(theta), sin(theta); -1*sin(theta), cos(theta)];
    y = R*v;
    y = y';
    y = [y, x(3)];
end

