function [y1] = dyn2(t,x1,F)
%State space dynamics of the tether for RK4 solver
    global G Mt M;
    %x1 = x';
   
    dist = norm([x1(1), x1(2), x1(3)],2);
%     L_cap = x/dist;
%     height = dist - R;
%     [lat, lon] = latlon(x(1:3));
%     F_b =  simpson_solver(Fb, L,L_cap, V, height, lat, lon, t, mu_r);
    
%     y = [x(4),x(5),x(6), -G*M*x(1)/(dist)^1.5 + F_b(1)/Mt ,...
%           -G*M*x(2)/(dist)^1.5 + F_b(2)/Mt , -G*M*x(3)/(dist)^1.5 + F_b(3)/Mt];
        
    %F = Fb(x1,t);
    y1 = [x1(4);x1(5);x1(6); ...
        -G*M*x1(1)/(dist)^3 + F(1)/Mt ;...
        -G*M*x1(2)/(dist)^3 + F(2)/Mt ; ...
        -G*M*x1(3)/(dist)^3 + F(3)/Mt];
    
    
end

