clc; close all; clear;

syms x x1 x2 x3 v v1 v2 v3 o o1 o2 o3 a a1 a2 a3 m u1 u2 u3 u4 l s g I1 I2 I3

x = [x1;x2;x3];
v = [v1;v2;v3];
o = [o1;o2;o3];
a = [a1;a2;a3];

I = [I1 0 0;
     0 I2 0;
     0 0 I3];

t_inv = [1 sin(a1)*tan(a2)  cos(a1)*tan(a2);
         0 cos(a1)          -sin(a1);
         0 sin(a1)/cos(a2)  cos(a1)/cos(a2)];

R_CE = [cos(a2)*cos(a3) sin(a1)*sin(a2)*cos(a3)-cos(a1)*sin(a3) sin(a1)*sin(a3)+cos(a1)*sin(a2)*cos(a3);
        cos(a2)*sin(a3) cos(a1)*cos(a3)+sin(a1)*sin(a2)*sin(a3) cos(a1)*sin(a2)*sin(a3)-sin(a1)*cos(a3);
        -sin(a2)        sin(a1)*cos(a2)                         cos(a1)*cos(a2)];

z = [x;a;v;o];
u = [u1; u2; u3; u4];

z0 = [x;zeros([9 1])];
u0 = [m*g/4; m*g/4; m*g/4; m*g/4];
z_dot = [v; t_inv*o; 1/m*R_CE*[0; 0; u1+u2+u3+u4]-[0;0;g]; I\([(u2 - u4)*l; (u3 - u1)*l; (u1 - u2 + u3 - u4)*s] - cross(o, I*o))];

df_dz = jacobian(z_dot, [x1; x2; x3; a1; a2; a3; v1; v2; v3; o1; o2; o3]);
df_dz = subs(df_dz, [x1;x2;x3;a1;a2;a3;v1;v2;v3;o1;o2;o3], [x1;x2;x3;0;0;0;0;0;0;0;0;0]);
df_dz = subs(df_dz, [u1; u2; u3; u4], [m*g/4; m*g/4; m*g/4; m*g/4]);
A = df_dz

df_du = jacobian(z_dot, [u1; u2; u3; u4;]);
df_du = subs(df_du, [x1;x2;x3;a1;a2;a3;v1;v2;v3;o1;o2;o3], [x1;x2;x3;0;0;0;0;0;0;0;0;0]);
df_du = subs(df_du, [u1; u2; u3; u4], [m*g/4; m*g/4; m*g/4; m*g/4]);
B = df_du

[I\[0; -l; s] I\[l; 0; -s] I\[0; l; s] I\[-l; 0; -s]];

z_lin = A*(z-z0) + B*(u - u0);

C = [B A*B A^2*B A^3*B A^4*B A^5*B A^6*B A^7*B A^8*B A^9*B A^10*B A^11*B];
rank(C)

% desired eigenvalues of -1, -2, -3, ... -11
B = subs(B, [g,l,m,I1,I2,I3,s], [9.81,0.2,0.5,1.24,1.24,2.48,0.01]);
A = subs(A, [g,l,m,I1,I2,I3,s], [9.81,0.2,0.5,1.24,1.24,2.48,0.01]);
B = double(B);
A = double(A);
desired_polynomial = [-1 -2 -3 -4 -5 -6 -7 -8 -9 -10 -11 -12];
% K = place(A, B, desired_polynomial);
Q = 1/60*eye(12);
R = eye(4);
[K,S,P] = lqr(A,B,Q,R);
K(K<0.001) = 0;


syms ex1 ex2 ex3 ev1 ev2 ev3 eo1 eo2 eo3 ea1 ea2 ea3 v_1 v_2 v_3 v_4 xd1 xd2 xd3
z_desired = [xd1; xd2; xd3;zeros([9 1])];
u_desired = u0;
e = z_desired - z;
v = u_desired - u;
e_dot = -subs(z_dot, z, z_desired - [ex1;ex2;ex3;ea1;ea2;ea3;ev1;ev2;ev3;eo1;eo2;eo3]);
e_dot = subs(e_dot, u, u_desired - [v_1; v_2; v_3; v_4]);

df_de = jacobian(e_dot, [ex1;ex2;ex3;ea1;ea2;ea3;ev1;ev2;ev3;eo1;eo2;eo3]);
df_de = subs(df_de, [ex1;ex2;ex3;ea1;ea2;ea3;ev1;ev2;ev3;eo1;eo2;eo3; v_1; v_2; v_3; v_4], zeros([16 1]));

df_dv = jacobian(e_dot, [v_1; v_2; v_3; v_4]);
df_dv = subs(df_dv, [ex1;ex2;ex3;ea1;ea2;ea3;ev1;ev2;ev3;eo1;eo2;eo3; v_1; v_2; v_3; v_4], zeros([16 1]));

e_lin = df_de*e + df_dv * v;

F = subs(df_dv, [g,l,m,I1,I2,I3,s], [9.81,0.2,0.5,1.24,1.24,2.48,0.01]);
E = subs(df_de, [g,l,m,I1,I2,I3,s], [9.81,0.2,0.5,1.24,1.24,2.48,0.01]);
F = double(F);
E = double(E);

Q = 1/30*[eye(6) zeros(6);
     zeros([6 12])];
R = eye(4);
[K2,S2,P2] = lqr(E,F,Q,R);
K2(K2<0.001) = 0;