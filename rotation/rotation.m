%%
function rotation()

function Rotx = Rx(data)
    Rotx = [1 0 0;0 cosd(data) -sind(data); 0 sind(data) cosd(data)];
end

function Roty = Ry(data)
    Roty = [cosd(data) 0 -sind(data); 0 1 0; sind(data) 0 cosd(data)];
end

function Rotz = Rz(data)
    Rotz = [cosd(data) -sind(data) 0; sind(data) cosd(data) 0; 0 0 1];
end

function qua = R2Q(R)
    q0 = sqrt(trace(R)+1)/2;
    q1 = (R(3,2)-R(2,3))/(4*q0);
    q2 = (R(1,3)-R(3,1))/(4*q0);
    q3 = (R(2,1)-R(1,2))/(4*q0);
    qua = [q0;q1;q2;q3];
end

function rot = Q2R(Q)
    q0=Q(1); q1=Q(2); q2=Q(3); q3=Q(4);
    rot = [1-2*q2*q2-2*q3*q3 2*q1*q2-2*q0*q3   2*q1*q3+2*q0*q2; ...
           2*q1*q2+2*q0*q3   1-2*q1*q1-2*q3*q3 2*q2*q3-2*q0*q1; ...
           2*q1*q3-2*q0*q2   2*q2*q3+2*q0*q1   1-2*q1*q1-2*q2*q2];
end

function qua = Qmulti(Q1,Q2)
    a1=Q1(1); b1=Q1(2); c1=Q1(3); d1=Q1(4);
    a2=Q2(1); b2=Q2(2); c2=Q2(3); d2=Q2(4);
    qua = [a1*a2-b1*b2-c1*c2-d1*d2;...
           a1*b2+b1*a2+c1*d2-d1*c2;...
           a1*c2-b1*d2+c1*a2+d1*b2;...
           a1*d2+b1*c2-c1*b2+d1*a2];
end

function qua = Qconj(Q)
    qua = [Q(1);-Q(2);-Q(3);-Q(4)];
end

function Ome = Omega(angle)
    x=angle(1); y=angle(2); z=angle(3);
    Ome = [0 -x -y -z;...
           x  0  z -y;...
           y -z  0  x;...
           z  y -x  0];
end

function em = Skew(angle)
    n=angle(1); e=angle(2); d=angle(3);
    em = [ 0 -d  e;...
           d  0 -n;...
          -e  n  0;];
end

%%

% Pw = Rwb1*Rb1b2*Pb2
Pw=[-2;0;0];
Pb1=[-1;sqrt(3);0];
Pb2=[1;sqrt(3);0];

Rwb1 = Rz(60) *Ry(0) * Rx(0);
Rb1b2 = Rz(60) *Ry(0) * Rx(0);

Rb1b2*Pb2; % == Pb1
Rwb1*Rb1b2*Pb2; % == Pb2

%%


Rwb1;
Rb1b2;
Rwb2=Rwb1*Rb1b2;

Qwb1=R2Q(Rwb1);
Qb1b2=R2Q(Rb1b2);

% Qwb1*Qb1b2 = R2Q(Rwb1*Rb1b2)
Qwb2=R2Q(Rwb2);
Qmulti(Qwb1,Qb1b2);

% Rwb2 = Q2R(Qwb1*Qb1b2)
Rwb2;
Q2R(Qmulti(Qwb1,Qb1b2));

%%

Qmulti(Qwb2,Qmulti([0;Pb2],Qconj(Qwb2)));

%%

angle=[0.001;0.002;0.0015];
Q2R((diag(ones(4,1))+0.5*Omega(angle))*Qwb2);
(diag(ones(3,1))-Skew(angle))*Rwb2;

%%

a = [3;4;5];
Skew(Rwb2*a)

Rwb2*Skew(a)*Rwb2'

end