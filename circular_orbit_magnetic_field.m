clear;
clc;

size=0.001 ;
i=0 ;
f=1000;
t = i : size : f;


vx = zeros(1,(f-i)/size +1,'double');
vy = zeros(1,(f-i)/size +1,'double');
vz = zeros(1,(f-i)/size +1,'double');
x=zeros(1,(f-i)/size +1,'double');
y=zeros(1,(f-i)/size +1,'double');
z=zeros(1,(f-i)/size +1,'double');

vx(1)=0;
vy(1)=7.3176e3;
vz(1)=0;

x(1)=5.1e5+6.4e6;                  % limit of igrfmagm is 600 km as height
y(1)=0;
z(1)=0;

G=6.67e-11;
M=5.7e24;          %kg        %initialization
Mt=10   %kg


%V=[vx vy vz];
R=6.4e6
r=[x y z];
L=100;
d=1

xo=[1 0 0];
yo=[0 1 0];
zo=[0 0 1];


%%


for d=1:(f-i)/size
V=[vx(d) vy(d) vz(d)];
dist=((x(d)^2+y(d)^2+z(d)^2)^0.5);
                                              % lat,long
height=dist-R;
if z(d)==0
    theta=0
else
theta=(z(d)/abs(z(d)))*(acos(((x(d)^2+y(d)^2)^0.5)/((x(d)^2+y(d)^2+z(d)^2)^0.5)))*90/(pi/2)   ; %lat
end
% alpha longit.
if y(d)==0
    if x(d)>=0 
        alpha = 0;
    else
        alpha = pi;
    end
    
     
else
  alpha=(y(d)/abs(y(d)))*acos(x(d)/((x(d)^2+y(d)^2)^0.5))*90/(pi/2)    ; %x,y,z>0 
end
%x axis is 0-longitude and latitutde


%[magFieldVector,horIntensity,declination,inclination,totalIntensity,magFieldSecVariation,secVariationHorizontal,secVariationDeclination,secVariationInclination,secVariationTotal] = igrfmagm(height,theta,alpha,decyear(2015,7,4),12);

[B, hi, de, in, ti, mf, svh, svd, svi, svt] = igrfmagm(height, theta, alpha, decyear(2016,5,1),12 );
L_vector=(L/((x(d)^2+y(d)^2+z(d)^2)^0.5))*[x(d) y(d) z(d)];
%%
%=magFieldVector;
                           %convert from geod to ecif
ur=cross(V,B)
ut=dot(L_vector,ur)


uk=cross(L_vector,B);   %magnetic field taken as a constant??????
    
   
        
ax=-G*M*x(d)/(x(d)^2+y(d)^2)^1.5 + (ut/(L*Mt))*dot(uk,xo) ;
ay=-G*M*y(d)/(x(d)^2+y(d)^2)^1.5 + (ut/(L*Mt))*dot(uk,yo) ;
az=-G*M*z(d)/(x(d)^2+y(d)^2)^1.5 + (ut/(L*Mt))*dot(uk,zo) ;

vx(d+1)=vx(d)+ ax*(size);
vy(d+1)=vy(d)+ ay*(size);
vz(d+1)=vz(d)+ az*(size);

x(d+1)= x(d) + vx(d)*(size);
y(d+1)= y(d) + vy(d)*(size);
z(d+1)= z(d) + vz(d)*(size);

end
for r=1:(f-i)/size +1
    m(r)=(x(r)^2 + y(r)^2 + z(r)^2)^0.5;
end

plot(t,m)
xlabel('time')

