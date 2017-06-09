%tether simulation assuming : circular orbit, tether radially pointing

clear;
clc;

step_size=0.001 ;
time_i=0 ;
time_f=1;
time = time_i : step_size : time_f;


vx = zeros(1,(time_f-time_i)/step_size +1,'double');
vy = zeros(1,(time_f-time_i)/step_size +1,'double');
vz = zeros(1,(time_f-time_i)/step_size +1,'double');
x=zeros(1,(time_f-time_i)/step_size +1,'double');
y=zeros(1,(time_f-time_i)/step_size +1,'double');
z=zeros(1,(time_f-time_i)/step_size +1,'double');

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
V=[vx(n) vy(n) vz(n)];
dist=((x(n)^2+y(n)^2+z(n)^2)^0.5);
                                              % lat,long
height=dist-R;
if z(n)==0
    theta=0
    else
    theta=(z(n)/abs(z(n)))*(acos(((x(n)^2+y(n)^2)^0.5)/((x(n)^2+y(n)^2+z(n)^2)^0.5)))*90/(pi/2)   ; %lat
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


%[magfieldVector,horIntensity,declination,inclination,totalIntensity,magfieldSecVariation,secVariationHorizontal,secVariationDeclination,secVariationInclination,secVariationTotal] = igrfmagm(height,theta,alpha,decyear(2015,7,4),12);

[B] = igrfmagm(height, theta, alpha, decyear(2016,5,1),12 );
L_vector=(L/((x(n)^2+y(n)^2+z(n)^2)^0.5))*[x(n) y(n) z(n)];
B = B*1e-9; %convert nanotesla to tesla
%=magfieldVector;
%convert from geod to ecitime_f


VxB=cross(V,B);
emf=dot(L_vector,VxB);


L_vectorxB=cross(L_vector,B);   %magnetic field taken as a constant
    
   
        
ax=-G*M*x(n)/(x(n)^2+y(n)^2 + z(n)^2)^1.5 + (emf/(L*Mt))*dot(L_vectorxB,xo) ;
ay=-G*M*y(n)/(x(n)^2+y(n)^2 + z(n)^2)^1.5 + (emf/(L*Mt))*dot(L_vectorxB,yo) ;
az=-G*M*z(n)/(x(n)^2+y(n)^2 + z(n)^2)^1.5 + (emf/(L*Mt))*dot(L_vectorxB,zo) ;

vx(n+1)=vx(n)+ ax*(step_size);
vy(n+1)=vy(n)+ ay*(step_size);
vz(n+1)=vz(n)+ az*(step_size);

x(n+1)= x(n) + vx(n)*(step_size);
y(n+1)= y(n) + vy(n)*(step_size);
z(n+1)= z(n) + vz(n)*(step_size);

end

%calculate m after every instant, where m is radial distance

for r=1:(time_f-time_i)/step_size +1
    radius(r)=(x(r)^2 + y(r)^2 + z(r)^2)^0.5;
end

plot(time,radius)
xlabel('time')

