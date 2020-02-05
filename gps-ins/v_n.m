syms n e d;
Rs(n,e,d) = [0 -d e; d 0 -n; -e n 0];
syms w_e B L h M N v_n v_e v_d;

w_ie = [w_e*cos(B); 0; -w_e*sin(B)];
w_en = [v_e/(N+h);-v_n/(M+h);-v_e*tan(B)/(N+h)];
v_n = [v_n; v_e; v_d];

R1=Rs(2*w_ie(1), 2*w_ie(2), 2*w_ie(3));
R2=Rs(w_en(1),w_en(2),w_en(3));

R3=(R1+R2)*v_n;

%%
clc;
diff(R3,B)
diff(R3,L)
diff(R3,h)
