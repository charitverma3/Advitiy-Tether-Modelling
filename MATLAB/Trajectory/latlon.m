function [ lat,lon ] = latlon( x)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    if x(3)==0
        lat=0;
    else
        lat=(x(3)/abs(x(3)))*(acos(((x(1)^2+x(2)^2)^0.5)/((x(1)^2+x(2)^2+x(3)^2)^0.5)))*90/(pi/2); 
    end
    
    % longitude calculation given position, lon is longitude
    if x(2)==0
        if x(1)>=0 
            lon = 0;
        else
            lon = pi;
        end
    else
      lon=(x(2)/abs(x(2)))*acos(x(1)/((x(1)^2+x(2)^2)^0.5))*90/(pi/2);    ; %x,y,z>0 
    end
    %x axis is intersection of 0 longitude and 0 latitutde

end

