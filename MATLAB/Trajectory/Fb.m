function [ F,e1] = Fb(x,t)
%Fb calculates magnetic force on tether, using Simpson's rule for
%integration, x is the state
%getting mangetic field from IGRF model
    global nL L R mu_r day E;
    
    dL = L/nL;
    pos = [x(1) x(2) x(3)];
    v_i = [x(4) x(5) x(6)];
    v_e = ecif2ecef(v_i,t);
    dist = norm(pos,2);
    dL_cap = -pos/dist;
    dL_cap_e = ecif2ecef(dL_cap,t);
    pos1 = ecif2ecef(pos,t);
    [lat, lon] = latlon(pos1);
    dL_vector = dL_cap*dL;
    dL_vector_e = ecif2ecef(dL_vector,t);
    F = [0,0,0];
    height = dist - R;
    e1 = 0;
    %height
    for i = 1:nL
        if nL == 1
            F = dF_b(height);
        else
            %F = F + (dF_b(height) + 4*dF_b(height + dL/2) + dF_b(height + dL))/6;
        end
        height = height - nL*dL;
    end
    
    
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


    
    function [dF_i] = dF_b(height)
        %tic
        B = igrf1(day, lat, lon, height/1e3,'geod'); %height in km
        
        %toc
        B1 = B*1e-9; %convert from nanotesla to tesla
        B = ned2ecef1(B1,lat,lon,height);
        e1 = dot(dL_cap_e, cross(v_e, B)); %thisfor is actually emf/L 
        dF_e = e1*(cross(dL_vector_e,B));
        dF_e = dF_e/mu_r;
        dF_i = ecef2ecif(dF_e,t);
%         e1 = dot(dL_cap, cross(v_i, B)); %this is actually emf/L 
%         dF_i = e1*(cross(dL_vector,B));
%         dF_i = dF_i/mu_r;
        
    end
    
end