%tether simulation assuming : circular orbit, tether radially pointing.
%used matrices instead of vectors for position and velocities. assumed variable magnetic field over
%tether, did euler integration for tether and used ode45 for orbit.

clear;
clc;

step_size=0.01; %time step for simulation
time_i=0; %start time of simulation
time_f=1; %end time of simulation
time = time_i : step_size : time_f; %array of time
nL = 5; %number of parts length of tether is divided for euler integration of force


%vx = zeros(1,(time_f-time_i)/step_size +1,'double');
%vy = zeros(1,(time_f-time_i)/step_size +1,'double');
%vz = zeros(1,(time_f-time_i)/step_size +1,'double');
x=zeros(1,(time_f-time_i)/step_size +1,'double'); %x co-ordinate of position w.r.t earth
y=zeros(1,(time_f-time_i)/step_size +1,'double'); %y co-ordinate of position w.r.t earth
z=zeros(1,(time_f-time_i)/step_size +1,'double'); %z co-ordinate of position w.r.t earth


V = zeros((time_f-time_i)/step_size +1, 3, 'double'); %velocity as a functoin of time
pos = zeros((time_f-time_i)/step_size +1, 3, 'double'); %position as a functoin of time
%V_cm = zeros((time_f-time_i)/step_size +1, 3, 'double');
%pos_cm = zeros((time_f-time_i)/step_size +1, 3, 'double');
%eps = zeros(1,(time_f-time_i)/step_size +1,'double');

x(1)=5.1e5+6.4e6;                  % limit of igrtime_fmagm is 600 km as height
y(1)=0;
z(1)=0;
pos(1,:) = [x(1), y(1), z(1)];

G=6.67e-11; %universal gravitational constant, SI
M=5.7e24; %mass of earth, kg
Mt=10;   %total mass of system, kg

vx(1)=0;
vy(1)=sqrt(G*M/x(1));
vz(1)=0;

V(1,:)=[vx(1) vy(1) vz(1)];
R=6.4e6; %radius of earth, m
r=[x y z];
L=100; %length of tether, m
n=1; 
dL = L/nL; 
rho = 0.05; %resistance per unit length of tether, ohm/m
xo=[1 0 0]; %unit vector in x direction
yo=[0 1 0]; %unit vector in y direction
zo=[0 0 1]; %unit vector in z direction

times_output_log=[]; %for RK4(ode45) solver
soln_output_log=zeros(6,(time_f-time_i)/step_size); %for RK4(ode45) solver
%%


array_a1=zeros(1,(time_f-time_i)/step_size);
array_a2=zeros(1,(time_f-time_i)/step_size);


%this part to find numerical error


for n=1:(time_f-time_i)/step_size

    
    %[x(n) y(n) z(n)] = pos(n,:);
    pos1 = [x(n) y(n) z(n)]; %defined for convenience
    dist = norm(pos1,2); 
    %dist_cm = norm(pos_cm(n,:),2);
    %eps(n) = dist_cm - norm([x(1) y(1) z(1)],2);
    height=dist-R; %height of satellite from earth's surface
    
    %latitude calculation given position, theta is latitude
    if z(n)==0
        theta=0;
        else
        theta=(z(n)/abs(z(n)))*(acos(((x(n)^2+y(n)^2)^0.5)/((x(n)^2+y(n)^2+z(n)^2)^0.5)))*90/(pi/2); 
    end
    
    % longitude calculation given position, alpha is longitude
    if y(n)==0
        if x(n)>=0 
            alpha = 0;
        else
            alpha = pi;
        end
    else
      alpha=(y(n)/abs(y(n)))*acos(x(n)/((x(n)^2+y(n)^2)^0.5))*90/(pi/2);    ; %x,y,z>0 
    end
    %x axis is intersection of 0 longitude and 0 latitutde

    L_vector=L*pos1/dist; %length vector of tether     

    l=0;
    F = 0;
    for i=1:nL
        dL_cap = L_vector/L;
        dL_vector = dL_cap*dL;
        B = igrfmagm(height - l, theta, alpha, decyear(2016,5,1) + n*step_size,12 );
        B = B*1e-9; %convert nanotesla to tesla
        %B = B*0; %inserted fot testing purpose
        dF = dot(dL_cap, cross(V(n,:), B))*(cross(dL_vector,B));
        dF = dF/rho;
        F = F + dF;
        l = l + nL*dL;
    end
    %solns is the vector of (x,y,z,xdot,ydot,zdot)and times is time array
    %for rk4
    %F = F*0;
    f=@(times,soln)[soln(4);soln(5);soln(6); -G*M*soln(1)/(soln(1)^2+soln(2)^2+soln(3)^2)^1.5 + dot(F,xo)/Mt ;-G*M*soln(2)/(soln(1)^2+soln(2)^2+soln(3)^2)^1.5 + dot(F,yo)/Mt ; -G*M*soln(3)/(soln(1)^2+soln(2)^2+soln(3)^2)^1.5 + dot(F,zo)/Mt];
    [times_output,soln_output]=ode45(f,[n*step_size, (n+1)*step_size],[x(n),y(n),z(n),vx(n),vy(n),vz(n)]);
    %initial conditions for every iteration given
    
    %creating log of times_ouput and soln_output
    %times_output_log(n)=times_output(n);
    %soln_output_log(:,n)=soln_output(n);
    
    x(n+1)=soln_output(end,1);
    y(n+1)=soln_output(end,2);
    z(n+1)=soln_output(end,3);
    vx(n+1)=soln_output(end,4);
    vy(n+1)=soln_output(end,5);
    vz(n+1)=soln_output(end,6);
    
   % a_cm = -G*M*pos_cm/dist_cm^3
   % V_cm(n+1,:) = V_cm(n,:) + a_cm*step_size;
    %pos_cm(n+1,:) = pos_cm(n,:) + V_cm(n,:)*step_size;
    
    %a = -G*M*pos1/dist^3 + F/Mt; 
    %V(n+1,:) = V(n,:) + a*step_size; 
    %pos(n+1,:) = pos(n,:) + V(n,:)*step_size;
    %x(n+1) = pos(n+1,1);
    %y(n+1) = pos(n+1,2);
    %z(n+1) = pos(n+1,3);
end

%calculate m after every instant, where m is radial distance
%%
for r=1:(time_f-time_i)/step_size +1
    radius(r)=(x(r)^2 + y(r)^2 + z(r)^2)^0.5;
end

plot(time,radius)
xlabel('time')





