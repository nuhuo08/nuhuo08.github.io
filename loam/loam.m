clc;clear;

%%
% Rotation matrix
syms data;
Rx(data) = [1 0 0;0 cos(data) -sin(data); 0 sin(data) cos(data)];
Ry(data) = [cos(data) 0 sin(data); 0 1 0; -sin(data) 0 cos(data)];
Rz(data) = [cos(data) -sin(data) 0; sin(data) cos(data) 0; 0 0 1];

%%

syms r o x;

Rnb=Rz(x)*Rx(o)*Ry(r);
Rbn=Ry(-r)*Rx(-o)*Rz(-x);

%%

syms r p y;

R=Rz(y)*Rx(p)*Ry(r);

q0=sqrt(1+R(1,1)+R(2,2)+R(3,3))/2
q1=(R(3,2)-R(2,3))/4/q0
q2=(R(1,3)-R(3,1))/4/q0
q3=(R(2,1)-R(1,2))/4/q0


%%

% Optimization equation
syms s;
syms rx ry rz;
% Rbw = Rz(-rz) * Rx (-rx) * Ry(-ry);
Rbw = Ry(-s*ry) * Rx(-s*rx) *  Rz(-s*rz);
% Rbw = Ry(ry) * Rx(rx) *  Rz(rz);

arx = diff(Rbw,rx)
ary = diff(Rbw,ry);
arz = diff(Rbw,rz);

%%
% Accumulate Rotaion
syms x1 y1 z1 x2 y2 z2;
Rwb1 = Ry(y1) * Rx(x1) *  Rz(z1);
Rbw2 = Ry(y2) * Rx(x2) *  Rz(z2);

Rwb3 = Rwb1*Rbw2

syms x y z;
Rwb3 = Ry(y) * Rx(x) *  Rz(z)

%%
% Plugin IMU
clc;
syms alx aly alz blx bly blz bcx bcy bcz;

R1 = Ry(alz) * Rx(alx) *  Rz(aly);
R2 = Rz(-bly) * Rx(-blx) * Ry(-blz);
% R2 = Ry(bly) * Rx(blx) *  Rz(blz);
R3 = Ry(bcz) * Rx(bcx) *  Rz(bcy);

R = R1 * R2 * R3;
pretty(R(2,3))

%%

% Point to Line
clc; clear;

syms x0 y0 z0 x1 y1 z1 x2 y2 z2;


a012 = sqrt(((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1)) ...
         * ((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1))  ...
         + ((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1)) ...
         * ((x0 - x1)*(z0 - z2) - (x0 - x2)*(z0 - z1))  ...
         + ((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1)) ...
         * ((y0 - y1)*(z0 - z2) - (y0 - y2)*(z0 - z1)));

l12 = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));

f = a012/l12;

diff(f,x0)

pretty(diff(f,x0))

%% TODO
% Optimize Axis-Angle

syms T4 T5 T6 x1 y1 z1 x2 y2 z2;

syms theta omega_1 omega_2 omega_3 omega_hat omega_hat2 R;

theta = sqrt(T4 * T4 + T5 * T5 + T6 * T6);

omega_1 = T4 / theta;
omega_2 = T5 / theta;
omega_3 = T6 / theta;

omega_hat = [0 -omega_3 omega_2; omega_3 0 -omega_1; -omega_2 omega_1 0];
omega_hat2 = omega_hat^2;

R = eye(3,3) + omega_hat * sin(theta) + omega_hat2*(1-cos(theta));

x2 = R(1,1) * x1 + R(1,2) * y1 + R(1,3) * z1;
y2 = R(2,1) * x1 + R(2,2) * y1 + R(2,3) * z1;
z2 = R(3,1) * x1 + R(3,2) * y1 + R(3,3) * z1;


