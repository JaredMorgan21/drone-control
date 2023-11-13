function dz = quadrotor(t, z, u, p, r, n, zd)
% State vector definition
%
%      x1, x2, x3, phi, theta, psi, dx1, dx2, dx3, omega1, omega2, omega3
% z = [z1, z2, z3,  z4,    z5,  z6,  z7,  z8,  z9,    z10,    z11,    z12]
%
% Parameter vector definition
%
%       g,  l,  m, I11, I22, I33, mu, sigma
% p = [p1, p2, p3,  p4,  p5,  p6, p7,    p8]


% Forming the moment of inertia tensor based on the parametr values
I = diag(p(4:6)); 


% Rotation matrix mapping body fixed frame C to inertial frame E
R = [ cos(z(5))*cos(z(6)), sin(z(4))*sin(z(5))*cos(z(6)) - cos(z(4))*sin(z(6)), sin(z(4))*sin(z(6)) + cos(z(4))*sin(z(5))*cos(z(6));
      cos(z(5))*sin(z(6)), cos(z(4))*cos(z(6)) + sin(z(4))*sin(z(5))*sin(z(6)), cos(z(4))*sin(z(5))*sin(z(6)) - sin(z(4))*cos(z(6));
               -sin(z(5)),                                 sin(z(4))*cos(z(5)),                                 cos(z(4))*cos(z(5))];

% f = [z(7);
%      z(8);
%      z(9);
%      z(10) + z(12)*cos(z(4))*tan(z(5)) + z(11)*sin(z(4))*tan(z(5));
%      z(11)*cos(z(4)) - z(12)*sin(z(4));
%      (z(12)*cos(z(4)))/cos(z(5)) + (z(11)*sin(z(4)))/cos(z(5));
%      0;
%      0;
%      -p(1);
%      -(p(5)*z(11)*z(12) - p(6)*z(11)*z(12))/p(4);
%      (p(4)*z(10)*z(12) - p(6)*z(10)*z(12))/p(5);
%      -(p(4)*z(10)*z(11) - p(5)*z(10)*z(11))/p(6)];
% g = [                                        0,                                              0,                                              0,                                              0;
%                                              0,                                              0,                                              0,                                              0;
%                                              0,                                              0,                                              0,                                              0;
%                                              0,                                              0,                                              0,                                              0;
%                                              0,                                              0,                                              0,                                              0;
%                                              0,                                              0,                                              0,                                              0;
%  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3),  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3),  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3),  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3);
% -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3), -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3), -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3), -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3);
%                            (cos(z(4))*cos(z(5)))/p(3),                            (cos(z(4))*cos(z(5)))/p(3),                            (cos(z(4))*cos(z(5)))/p(3),                            (cos(z(4))*cos(z(5)))/p(3);
%                                              0,                                           p(2)/p(4),                                              0,                                          -p(2)/p(4);
%                                          -p(2)/p(5),                                              0,                                           p(2)/p(5),                                              0;
%                                           p(8)/p(6),                                          -p(8)/p(6),                                           p(8)/p(6),                                          -p(8)/p(6)];
% 
% A =     [0         0         0         0         0         0    1.0000         0         0         0         0         0;
%          0         0         0         0         0         0         0    1.0000         0         0         0         0;
%          0         0         0         0         0         0         0         0    1.0000         0         0         0;
%          0         0         0         0         0         0         0         0         0    1.0000         0         0;
%          0         0         0         0         0         0         0         0         0         0    1.0000         0;
%          0         0         0         0         0         0         0         0         0         0         0    1.0000;
%          0         0         0         0    9.8100         0         0         0         0         0         0         0;
%          0         0         0   -9.8100         0         0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0         0         0         0         0         0;
%          0         0         0         0         0         0         0         0         0         0         0         0];
% 
% B =     [0         0         0         0;
%          0         0         0         0;
%          0         0         0         0;
%          0         0         0         0;
%          0         0         0         0;
%          0         0         0         0;
%          0         0         0         0;
%          0         0         0         0;
%     2.0000    2.0000    2.0000    2.0000;
%          0    0.1613         0   -0.1613;
%    -0.1613         0    0.1613         0;
%     0.0040   -0.0040    0.0040   -0.0040];
% 
% K = [0	0	2.73860000000000	0	0	3.16230000000000	0	0	0.966800000000000	0	0	19.8084000000000;
% 0	0	2.73860000000000	29.3486000000000	0	0	0	0	0.966800000000000	13.5078000000000	0	0;
% 2.23610000000000	0	2.73860000000000	0	29.3486000000000	3.16230000000000	3.72550000000000	0	0.966800000000000	0	13.5078000000000	19.8084000000000;
% 0	2.23610000000000	2.73860000000000	0	0	0	0	3.72550000000000	0.966800000000000	0	0	0];
% 
% 
% w = A*z + B*-K*(z-zd);
% % v1=1;v2=1;v3=1;
% % w = [0;
% %      0;
% %      0;
% %      zeros(9,1)];
% % u = -K*(z - zd) + p(3) * p(1) / 4
% u = g\(w - f)

% [u;z(5)]

% Adjusting thrust output based on feasible limits
u = max(min(u, p(7)), 0);

% Computing temporrary variables

% rt = torque vector induced by rotor thrusts
rt = [                   ( u(2) - u(4) )*p(2); 
                         ( u(3) - u(1) )*p(2); 
           ( u(1) - u(2) + u(3) - u(4) )*p(8)];

% Computing time derivative of the state vector
dz(1:3,1) = z(7:9,1);

dz(4:6,1) = [ z(10) + z(12)*cos(z(4))*tan(z(5)) + z(11)*sin(z(4))*tan(z(5));
                                          z(11)*cos(z(4)) - z(12)*sin(z(4));
                              (z(12)*cos(z(4)) + z(11)*sin(z(4)))/cos(z(5))];
                      
dz(7:9,1) = R*([0; 0; sum(u)] + r)/p(3) - [0; 0; p(1)];

dz(10:12,1) = I\( rt + n - cross( z(10:12,1) , I * z(10:12,1) ) );