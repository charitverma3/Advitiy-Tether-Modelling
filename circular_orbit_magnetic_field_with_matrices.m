%tether simulation assuming : circular orbit, tether radially pointing.
%used matrices instead of vectors

clear;
clc;

step_size=0.01 ;
time_i=0 ;
time_f=1;
time = time_i : step_size : time_f;
nL = 100;
dL = L/nL;

%vx = zeros(1,(time_f-time_i)/step_size +1,'double');
%vy = zeros(1,(time_f-time_i)/step_size +1,'double');
%vz = zeros(1,(time_f-time_i)/step_size +1,'double');
x=zeros(1,(time_f-time_i)/step_size +1,'double');
y=zeros(1,(time_f-time_i)/step_size +1,'double');
z=zeros(1,(time_f-time_i)/step_size +1,'double');

V = zeros((time_f-time_i)/step_size +1, 3, 'double');
pos = zeros((time_f-time_i)/step_size +1, 3, 'double');

vx(1)=0;
vy(1)=7.3176e3;
vz(1)=0;

x(1)=5.1e5+6.4e6;                  % limit of igrtime_fmagm is 600 km as height
y(1)=0;
z(1)=0;

G=6.67e-11;
M=5.7e24;          %kg        %initialization
Mt=10   %kg


%V=[vx vy vz];
R=6.4e6
r=[x y z];
L=100;
n=1

xo=[1 0 0];
yo=[0 1 0];
zo=[0 0 1];


%%


for n=1:(time_f-time_i)/step_size

    l=0;
    F = 0;
    [x(n) y(n) z(n)] = pos(n,:);
    pos1 = [x(n) y(n) z(n)];
    dist = norm(pos1,2);
                                                  
    height=dist-R;
    if z(n)==0
        theta=0
        else
        theta=(z(n)/abs(z(n)))*(acos(((x(n)^2+y(n)^2)^0.5)/((x(n)^2+y(n)^2+z(n)^2)^0.5)))*90/(pi/2); %lat
    end
    % alpha longit.
    if y(n)==0
        if x(n)>=0 
            alpha = 0;
        else
            alpha = pi;
        end
    else
      alpha=(y(n)/abs(y(n)))*acos(x(n)/((x(n)^2+y(n)^2)^0.5))*90/(pi/2)    ; %x,y,z>0 
    end
    %x axis is 0-longitude and latitutde

    L_vector=L*pos1/dist;
     %convert nanotesla to tesla



    for i=1:nL
        dL_cap = L_vector/L;
        dL_vector = dL_cap*dL
        B = igrfmagm(height - l, theta, alpha, decyear(2016,5,1),12 );
        B = B*1e-9; %convert nanotesla to tesla
        dF = dot(dL_cap, cross(V, B))*(cross(dL_vector,B));
        dF = dF/rho;
        F = F + dF;
        l = l + nL*dL;
    end


    a = -G*M*pos1/dist^3 + F/M; 
    V(n+1,:) = V(n,:) + a*time_step; 
    pos(n+1,:) = pos(n,:) + V(n+1,:)*time_step;
     

end

%calculate m after every instant, where m is radial distance

for r=1:(time_f-time_i)/step_size +1
    radius(r)=(x(r)^2 + y(r)^2 + z(r)^2)^0.5;
end

plot(time,radius)
xlabel('time')





