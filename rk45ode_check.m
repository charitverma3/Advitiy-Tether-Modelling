%code to try ode45 solver on equations of circular motion

GM=40.02*1e+13;
R=6400000+600000;
V=sqrt(GM/R);
N = 1000000;
x = zeros(1,N+1);
y = zeros(1,N+1);
z = zeros(1,N+1);
vx = zeros(1,N+1);
vy = zeros(1,N+1);
vz = zeros(1,N+1);
x(1) = R; y(1) = 0; vx(1) = 0; vy(1) = V; z(1)=0; vz(1)=0;
for i=1:N
    f=@(t,a)[a(4);a(5);a(6);-GM*a(1)/(a(1)^2+a(2)^2+a(3)^2)^1.5;-GM*a(2)/(a(1)^2+a(2)^2+a(3)^2)^1.5;-GM*a(3)/(a(1)^2+a(2)^2+a(3)^2)^1.5];
    [t,n]=ode45(f,[(i-1)/N,(i)/N],[x(i),y(i),z(i),vx(i),vy(i),vz(i)]);
    x(i+1) = n(end,1);
    y(i+1) = n(end,2);
    vx(i+1) = n(end,4);
    vy(i+1) = n(end,5);
    z(i+1) = n(end, 3);
    vz(i+1) = n(end,6);
end
xf=n(end,1);
yf=n(end,2);
zf = n(end,3);
xyz = sqrt(xf^2 + yf^2 + zf^2);
R - xyz