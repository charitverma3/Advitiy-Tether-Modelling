function [y] = propagate(x,h,t)
%State space dynamics of the tether for RK4 solver
    global G Mt M;
    dist = norm([x(1), x(2), x(3)],2);
%     L_cap = x/dist;
%     height = dist - R;
%     [lat, lon] = latlon(x(1:3));
%     F_b =  simpson_solver(Fb, L,L_cap, V, height, lat, lon, t, mu_r);
    
%     y = [x(4),x(5),x(6), -G*M*x(1)/(dist)^1.5 + F_b(1)/Mt ,...
%           -G*M*x(2)/(dist)^1.5 + F_b(2)/Mt , -G*M*x(3)/(dist)^1.5 + F_b(3)/Mt];
    k1 = h*dyn(x);
    k2 = h*dyn(x+k1/2); 
    k3 = h*dyn(x+k2/2);
    k4 = h*dyn(x+k3);
    y = x + (k1 + 2*k2 + 2*k3 + k4)/6;
    
    function[y1]= dyn(x1)
        F = Fb(x1,t);
        y1 = [x1(4),x1(5),x1(6), ...
            -G*M*x1(1)/(dist)^3 + F(1)/Mt ,...
            -G*M*x1(2)/(dist)^3 + F(2)/Mt , ...
            -G*M*x1(3)/(dist)^3 + F(3)/Mt];
    end
end

