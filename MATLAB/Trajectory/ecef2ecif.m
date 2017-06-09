function [ y ] = ecef2ecif( x,t )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    v = [x(1) x(2)];
    v = v';
    global w_earth;
    theta = w_earth*t;
    R = [cos(theta), -1*sin(theta); sin(theta), cos(theta)];
    y = R*v;
    y = y';
    y = [y,x(3)];
end

