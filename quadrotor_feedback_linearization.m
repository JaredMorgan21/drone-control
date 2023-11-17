function u = quadrotor_feedback_linearization(tdiff, z, p, v)
%      x1, x2, x3, phi, theta, psi, dx1, dx2, dx3, omega1, omega2, omega3
% z = [z1, z2, z3,  z4,    z5,  z6,  z7,  z8,  z9,    z10,    z11,    z12]
%       g,  l,  m, I11, I22, I33, mu, sigma
% p = [p1, p2, p3,  p4,  p5,  p6, p7,    p8]

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
    % 
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
    % k = 1;
    % b = 1;
    % w = [zeros([6,1]); 0;0;-sin(t) + b*(cos(t) - z(9) + k*(sin(t)-z(3)));0;0;0];
    % % w = [zeros([6,1]);0;0;0;0;0;0];
    % 
    % u = g\(w - f);
    persistent xi;
    persistent zeta;

    if isempty(xi)
        xi = 0;
        zeta = 0;
    end

    v = [p(3)*p(1);0;0;0];

    delta = [(sin(z(4))*sin(z(6)) + cos(z(4))*sin(z(5))*cos(z(6)))/p(3) 0 0 0;
             (sin(z(4))*sin(z(6)) + cos(z(4))*sin(z(5))*cos(z(6)))/p(3) 0 0 0;
             (sin(z(4))*sin(z(6)) + cos(z(4))*sin(z(5))*cos(z(6)))/p(3) 0 0 0;
             0 0 sin(z(4))/(cos(z(5))*p(5)) cos(z(4))/(cos(z(5))*p(6))]
    b = [z(7);
         z(8);
         z(9);
         sin(z(4))/cos(z(5))*z(11)+cos(z(4))/cos(z(5))*z(12)];
    alpha = -delta\b;
    beta = delta;

    u_hat = alpha + beta\v;
    xi = xi + u_hat(1)*tdiff;
    zeta = zeta + xi*tdiff;
    
    u = [zeta;
         u_hat(2);
         u_hat(3);
         u_hat(4);];
end