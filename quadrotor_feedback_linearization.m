function u = quadrotor_feedback_linearization(t, z, p, zd, A, B, K)
    f = [z(7);
         z(8);
         z(9);
         z(10) + z(12)*cos(z(4))*tan(z(5)) + z(11)*sin(z(4))*tan(z(5));
         z(11)*cos(z(4)) - z(12)*sin(z(4));
         (z(12)*cos(z(4)))/cos(z(5)) + (z(11)*sin(z(4)))/cos(z(5));
         0;
         0;
         -p(1);
         -(p(5)*z(11)*z(12) - p(6)*z(11)*z(12))/p(4);
         (p(4)*z(10)*z(12) - p(6)*z(10)*z(12))/p(5);
         -(p(4)*z(10)*z(11) - p(5)*z(10)*z(11))/p(6)];
    
    g = [                                        0,                                              0,                                              0,                                              0;
                                                 0,                                              0,                                              0,                                              0;
                                                 0,                                              0,                                              0,                                              0;
                                                 0,                                              0,                                              0,                                              0;
                                                 0,                                              0,                                              0,                                              0;
                                                 0,                                              0,                                              0,                                              0;
     (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3),  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3),  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3),  (sin(z(4))*sin(z(6)) + cos(z(4))*cos(z(6))*sin(z(5)))/p(3);
    -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3), -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3), -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3), -(cos(z(6))*sin(z(4)) - cos(z(4))*sin(z(5))*sin(z(6)))/p(3);
                               (cos(z(4))*cos(z(5)))/p(3),                            (cos(z(4))*cos(z(5)))/p(3),                            (cos(z(4))*cos(z(5)))/p(3),                            (cos(z(4))*cos(z(5)))/p(3);
                                                 0,                                           p(2)/p(4),                                              0,                                          -p(2)/p(4);
                                             -p(2)/p(5),                                              0,                                           p(2)/p(5),                                              0;
                                              p(8)/p(6),                                          -p(8)/p(6),                                           p(8)/p(6),                                          -p(8)/p(6)];
    
    % k = 1;
    % b = 1;
    % w = [zeros([6,1]); 0;0;-sin(t) + b*(cos(t) - z(9) + k*(sin(t)-z(3)));0;0;0];
    w = [zeros([6,1]); 0;0;0;0;0;1];
    u = g\(w - f);
end