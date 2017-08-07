%tether simulation assuming : circular orbit, tether radially pointing.
%used matrices instead of vectors for position and velocities. assumed variable magnetic field over
%tether, did euler integration for tether and used ode45 for orbit.
%x axis is intersection of 0 degree latitude and logitude, z axis towards
%north pole

clear;
clc;
%%
%global variables
global step_size L nL nT G M R Mt mu_r w_earth day E;
step_size=0.1; %time step for simulation
time_i=0; %start time of simulation
time_f=86400*2; %end time of simulation
time = time_i : step_size : time_f; %array of time
nL = 1; %number of parts length of tether is divided for euler integration of force
nT = (time_f-time_i)/step_size;
w_earth = 7.2921159e-5;
G=6.67e-11; %universal gravitational constant, SI
M=5.972e24; %mass of earth, kg
R=6371.8e3; %radius of earth, m
Mt=10;   %total mass of system, kg
L=100; %length of tether, m
%n=1; 
dL = L/nL; 
mu_r = 0.05*2; %resistance per unit length of tether, ohm/m
%day = decyear(2017,5,30); %for igrfmagm
day = datenum(2017,5,30); %for igrf1
%E = referenceSphere('earth', 'm'); %for NED to ECEF
%%
F = [0,0,0];
%vx = zeros(1,(time_f-time_i)/step_size +1,'double');
%vy = zeros(1,(time_f-time_i)/step_size +1,'double');
%vz = zeros(1,(time_f-time_i)/step_size +1,'double');
x=zeros(1,nT +1,'double'); %x co-ordinate of position w.r.t earth
y=zeros(1,nT +1,'double'); %y co-ordinate of position w.r.t earth
z=zeros(1,nT +1,'double'); %z co-ordinate of position w.r.t earth


V = zeros(nT +1, 3, 'double'); %velocity as a function of time
pos = zeros(nT +1, 3, 'double'); %position as a function of time
state = zeros(nT +1, 6, 'double'); %[position, velocity] as a function of time
state_cm = zeros(nT +1, 6, 'double');
pos_cm = zeros(nT +1, 3, 'double');
r = zeros(nT +1, 1, 'double');
r_cm = zeros(nT +1, 1, 'double');
energy = zeros(nT +1, 1, 'double');
energy_cm = zeros(nT +1, 1, 'double');
p_dedt = zeros(nT, 1, 'double');
p_fv1 = zeros(nT, 1, 'double');
p_fv2 = zeros(nT, 1, 'double');
emf = zeros(nT, 1, 'double');
current = zeros(nT, 1, 'double');
%eps = zeros(1,(time_f-time_i)/step_size +1,'double');

x(1)=5.1e5+R;                  % limit of igrtime_fmagm is 600 km as height
y(1)=0;
z(1)=0;
pos(1,:) = [x(1), y(1), z(1)];
%velocities for 98 degree orbit
vx(1)=0;
vy(1)=sqrt(G*M/norm([x(1),y(1),z(1)],2))*cos(98*pi/180);
vz(1)=sqrt(G*M/norm([x(1),y(1),z(1)],2))*sin(98*pi/180);
%velocities for equatorial
% vx(1)=0;
% vy(1)=sqrt(G*M/norm([x(1),y(1),z(1)],2));
% vz(1)=0;
V(1,:)=[vx(1) vy(1) vz(1)];
state(1,:) = [x(1), y(1), z(1), vx(1), vy(1), vz(1)];
state_cm(1,:) = [x(1), y(1), z(1), vx(1), vy(1), vz(1)];
energy(1) = 0.5*Mt*(norm(state(1,4:6),2))^2 - G*M*Mt/norm(state(1,1:3),2);

xo=[1 0 0]; %unit vector in x direction
yo=[0 1 0]; %unit vector in y direction
zo=[0 0 1]; %unit vector in z direction

times_output_log=[]; %for RK4(ode45) solver
soln_output_log=zeros(6,(time_f-time_i)/step_size); %for RK4(ode45) solver
%%
tic
for n=1:nT
    if mod(n,1000)==0
        n/10
    end
    
    t = time(n);
    [F,emf(n)] = Fb(state(n,:),t);
    emf(n) = emf(n)*L;
    [~,soln] = ode45(@(t,soln)dyn2(t,soln,F), [n*step_size, (n+1)*step_size], state(n,:)');
    state(n+1,:) = soln(end,:);
    pos(n+1,:) = state(n+1,1:3);
    energy(n+1) = 0.5*Mt*(norm(state(n+1,4:6),2))^2 - G*M*Mt/norm(state(n+1,1:3),2);
    [~,soln1] = ode45(@(t,soln1)dyn2(t,soln1,[0,0,0]), [n*step_size, (n+1)*step_size], state_cm(n,:)');
    state_cm(n+1,:) = soln1(end,:);
    pos_cm(n+1,:) = state_cm(n+1,1:3);
    F_g = Fg(state(n,1:3));
    p_fv1(n) = dot(state(n,4:6),F+F_g);
    p_dedt(n) = (energy(n+1) - energy(n))/step_size;
    
    
end
toc
%%
for n=1:nT+1
    r_cm(n) = norm(state_cm(n,1:3),2);
    r(n) = norm(pos(n,:),2);
end

save 98degree2.mat
%%
% 
% 
% array_a1=zeros(1,(time_f-time_i)/step_size);
% array_a2=zeros(1,(time_f-time_i)/step_size);
% 
% 
% %this part to find numerical error
% 
% 
% for n=1:(time_f-time_i)/step_size
% 
%     
%     %[x(n) y(n) z(n)] = pos(n,:);
%     pos1 = [x(n) y(n) z(n)]; %defined for convenience
%     dist = norm(pos1,2); 
%     %dist_cm = norm(pos_cm(n,:),2);
%     %eps(n) = dist_cm - norm([x(1) y(1) z(1)],2);
%     height=dist-R; %height of satellite from earth's surface
%     
%     %latitude calculation given position, lat is latitude
%     if z(n)==0
%         lat=0;
%     else
%         lat=(z(n)/abs(z(n)))*(acos(((x(n)^2+y(n)^2)^0.5)/((x(n)^2+y(n)^2+z(n)^2)^0.5)))*90/(pi/2); 
%     end
%     
%     % longitude calculation given position, lon is longitude
%     if y(n)==0
%         if x(n)>=0 
%             lon = 0;
%         else
%             lon = pi;
%         end
%     else
%       lon=(y(n)/abs(y(n)))*acos(x(n)/((x(n)^2+y(n)^2)^0.5))*90/(pi/2);    ; %x,y,z>0 
%     end
%     %x axis is intersection of 0 longitude and 0 latitutde
% 
%     L_vector=L*pos1/dist; %length vector of tether     
% 
%     l=0;
%     F = 0;
%     for i=1:nL
%         dL_cap = L_vector/L;
%         dL_vector = dL_cap*dL;
%         B = igrfmagm(height - l, lat, lon, decyear(2016,5,1) + n*step_size,12 );
%         B = B*1e-9; %convert nanotesla to tesla
%         %B = B*0; %inserted fot testing purpose
%         dF = dot(dL_cap, cross(V(n,:), B))*(cross(dL_vector,B));
%         dF = dF/mu_r;
%         F = F + dF;
%         l = l + nL*dL;
%     end
%     %solns is the vector of (x,y,z,xdot,ydot,zdot)and times is time array
%     %for rk4
%     %F = F*0;
%     f=@(times,soln)[soln(4);soln(5);soln(6); -G*M*soln(1)/(soln(1)^2+soln(2)^2+soln(3)^2)^1.5 + dot(F,xo)/Mt ;-G*M*soln(2)/(soln(1)^2+soln(2)^2+soln(3)^2)^1.5 + dot(F,yo)/Mt ; -G*M*soln(3)/(soln(1)^2+soln(2)^2+soln(3)^2)^1.5 + dot(F,zo)/Mt];
%     [times_output,soln_output]=ode45(f,[n*step_size, (n+1)*step_size],[x(n),y(n),z(n),vx(n),vy(n),vz(n)]);
%     %initial conditions for every iteration given
%     
%     %creating log of times_ouput and soln_output
%     %times_output_log(n)=times_output(n);
%     %soln_output_log(:,n)=soln_output(n);
%     
%     x(n+1)=soln_output(end,1);
%     y(n+1)=soln_output(end,2);
%     z(n+1)=soln_output(end,3);
%     vx(n+1)=soln_output(end,4);
%     vy(n+1)=soln_output(end,5);
%     vz(n+1)=soln_output(end,6);
%     
%    % a_cm = -G*M*pos_cm/dist_cm^3
%    % V_cm(n+1,:) = V_cm(n,:) + a_cm*step_size;
%     %pos_cm(n+1,:) = pos_cm(n,:) + V_cm(n,:)*step_size;
%     
%     %a = -G*M*pos1/dist^3 + F/Mt; 
%     %V(n+1,:) = V(n,:) + a*step_size; 
%     %pos(n+1,:) = pos(n,:) + V(n,:)*step_size;
%     %x(n+1) = pos(n+1,1);
%     %y(n+1) = pos(n+1,2);
%     %z(n+1) = pos(n+1,3);
% end
% 
% %calculate m after every instant, where m is radial distance
% %%
% for r=1:(time_f-time_i)/step_size +1
%     radius(r)=(x(r)^2 + y(r)^2 + z(r)^2)^0.5;
% end
% 
% plot(time,radius)
% xlabel('time')





