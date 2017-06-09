function [ F ] = Fg(x1)
%Fb calculates magnetic force on tether, using Simpson's rule for
%integration, x is the state
%getting mangetic field from IGRF model
    global G M; 
    dist = norm(x1,2);
    F = [-G*M*x1(1)/(dist)^3 ,...
        -G*M*x1(2)/(dist)^3 , ...
        -G*M*x1(3)/(dist)^3 ];
    