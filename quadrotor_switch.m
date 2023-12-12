function dz = quadrotor_switch(t, z, u1,u2,p, r1,r2, n1,n2, zd)
% State vector definition
%
%      x1, x2, x3, phi, theta, psi, dx1, dx2, dx3, omega1, omega2, omega3
% z = [z1, z2, z3,  z4,    z5,  z6,  z7,  z8,  z9,    z10,    z11,    z12]
%
% Parameter vector definition
%
%       g,  l,  m, I11, I22, I33, mu, sigma
% p = [p1, p2, p3,  p4,  p5,  p6, p7,    p8]

thres_dist_capture = 0.2/2;

% Forming the moment of inertia tensor based on the parametr values
I = diag(p(4:6)); 

% Rotation matrix mapping body fixed frame C to inertial frame E
R = [ cos(z(5))*cos(z(6)), sin(z(4))*sin(z(5))*cos(z(6)) - cos(z(4))*sin(z(6)), sin(z(4))*sin(z(6)) + cos(z(4))*sin(z(5))*cos(z(6));
      cos(z(5))*sin(z(6)), cos(z(4))*cos(z(6)) + sin(z(4))*sin(z(5))*sin(z(6)), cos(z(4))*sin(z(5))*sin(z(6)) - sin(z(4))*cos(z(6));
               -sin(z(5)),                                 sin(z(4))*cos(z(5)),                                 cos(z(4))*cos(z(5))];

% Adjusting thrust output based on feasible limits

dist_capture = norm(z(1:3) - zd(1:3));
if dist_capture <= thres_dist_capture % if captured
    u = u2;
    r = r2; % Note: Disturbance should be randomly generated (20231212: Ryo)
    n = n2;
    flag_mode = 1; % captured
    %%%%%%%%%%%%%
    % Ryo
    %%%%%%
    % target position should be the origin. Generate trajectory?
    %%%%%%%%%%%%%
else
    u = u1;
    r = r1;
    n = n1;
    flag_mode = 2; % not captured
end

if flag_mode == 1
    str_mode = 'Captured';
elseif flag_mode == 2
    str_mode = 'Not captured';
end

disp('-------------')
str_disp = [str_mode,', Distance: ',sprintf('%.3f',dist_capture)];
disp(str_disp);

u = max(min(u, p(7)), 0);
disp(u);
disp('-------------')

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

% integrator term eta (Ryo: Is this term accumulated over time?)
dz(13:24,1) = z(1:12) - zd(1:12);