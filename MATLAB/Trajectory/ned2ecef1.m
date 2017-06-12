function [ y ] = ned2ecef1( x,lat,lon,height )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    v = [-x(3), -x(1), x(2)]; %convert to spherical polar r theta phi
    v = v';
    

    theta = -lat + 90; %in degree, polar angle
    phi = lon; %in degree, azimuthal angle
    

    DCM = [sind(theta)*cosd(phi), cosd(theta)*cosd(phi), -sind(phi);
    		sind(theta)*sind(phi), cosd(theta)*sind(phi), cosd(phi);
    		cosd(theta), -sind(theta), 0]; 
    %for spherical to cartesian

    y = DCM*v;
    y = y';
    
end

