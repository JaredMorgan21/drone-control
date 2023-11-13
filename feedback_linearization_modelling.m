clc;clear;close all
syms x x1 x2 x3 v v1 v2 v3 o o1 o2 o3 a a1 a2 a3 u1 u2 u3 u4 
% syms z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9) z(10) z(11) z(12)
syms g l m I1 I2 I3 s 
% syms p(1) p(2) p(3) p(4) p(5) p(6) p(8)
x = [x1;x2;x3];
a = [a1;a2;a3];
v = [v1;v2;v3];
o = [o1;o2;o3];

t_inv = [1 sin(a1)*tan(a2)  cos(a1)*tan(a2);
         0 cos(a1)          -sin(a1);
         0 sin(a1)/cos(a2)  cos(a1)/cos(a2)];
I = [I1 0 0;
     0 I2 0;
     0 0 I3];

f = [v;t_inv*o;[0; 0; -g]; I\cross(o, I*o)];

b = [zeros([6 4]);
    (sin(a1)*sin(a3)+cos(a1)*sin(a2)*cos(a3))/m * ones([1 4]);
    (cos(a1)*sin(a2)*sin(a3)-sin(a1)*cos(a3))/m * ones([1 4]);
    (cos(a1)*cos(a2))/m * ones([1 4]);
    0 l/I1 0 -l/I1;
    -l/I2 0 l/I2 0;
    s/I3 -s/I3 s/I3 -s/I3]

% f = subs(f, [g,l,m,I1,I2,I3,s], [p(1) p(2) p(3) p(4) p(5) p(6) p(8)]);
% b = subs(b, [g,l,m,I1,I2,I3,s], [p(1) p(2) p(3) p(4) p(5) p(6) p(8)]);
% 
% f = subs(f, [x1 x2 x3 v1 v2 v3 o1 o2 o3 a1 a2 a3], [z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9) z(10) z(11) z(12)]);
% b = subs(b, [x1 x2 x3 v1 v2 v3 o1 o2 o3 a1 a2 a3], [z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9) z(10) z(11) z(12)]);