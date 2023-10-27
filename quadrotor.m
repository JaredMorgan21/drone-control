function dz = quadrotor(t, z, p, r, n, zd)
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

K = [0	0	0.0645497224367903	0	0	0.0645497224367898	0	0	0.142492446381779	0	0	2.82989968882797;
     0	0	0.0645497224367903	5.98384078259183	0	0	0	0	0.142492446381779	6.09164560569660	0	0;
     0.0912870929175277	0	0.0645497224367903	0	5.98384078259184	0.0645497224367894	0.345975089262573	0	0.142492446381779	0	6.09164560569661	2.82989968882795;
     0	0.0912870929175276	0.0645497224367903	0	0	0	0	0.345975089262572	0.142492446381779	0	0	0];
u = -K*(z - zd) + p(3) * p(1) / 4

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