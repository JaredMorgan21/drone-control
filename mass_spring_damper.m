clc; clear; close all;
syms l1 l2
l1=1;l2=1;
h = (l1+l2)/3;
w = (l1 + l2)/2;
delta = (l1+l2)/10;

syms m1 m2 q1 q2 q_dot1 q_dot2 g t
m1 = 2;m2=1;g=9.81;
M = [(m1 + m2)*l1^2+m2*l2^2+2*l1*l2*m2*cos(q2)   m2*l2*(l1*cos(q2)+l2);
        m2*l2*(l1*cos(q2)+l2)                       m2*l2^2];
phi = [(m1+m2)*g*l1*cos(q1)+m2*l2*(g*cos(q1+q2)-l1*q_dot2*(2*q_dot1+q_dot2)*sin(q2));
        m2*l2*(l1*q_dot1^2*sin(q2)+g*cos(q1+q2))];

L = [l1*cos(q1)+l2*cos(q1+q2);
     l1*sin(q1)+l2*sin(q1+q2)];

J = jacobian(L, [q1,q2]);
J_dot = [-cos(q1)*q_dot1-cos(q1+q2)*(q_dot1+q_dot2)     -cos(q1+q2)*(q_dot1+q_dot2);
         -sin(q1)*q_dot1-sin(q1+q2)*(q_dot1+q_dot2)     -sin(q1+q2)*(q_dot1+q_dot2)];

a = linspace(-2,2,200);
wall = [delta*a;delta*sin(a*20)+w]';
k_wall = 10000;

N = 0;
q_dot = [q_dot1;q_dot2];
z_dot = -J*q_dot;
z_ddot = [N-z_dot(1)+10/(l1+l2)*(L(1)-11*(l1+l2)/10);
            -(l1+l2)/3*sin(t)];
v = J\(z_ddot-J_dot*q_dot)
u = M*v+phi;


tspan = linspace(0,10,200);
q0 = [pi/2;-pi/2;0;0];

[t,q] = ode45(@(t, q) dqdt(t, q), tspan, q0);
figure(1)
plot(t,q(:,1:2));
legend("q1", "q2")
figure(2)
plot(t,[cos(q(:,1))+cos(q(:,1)+q(:,2)) sin(q(:,1))+sin(q(:,1)+q(:,2))]);
legend("x", "y")
figure(3)

plot(wall(2,:), wall(1,:))
for i = 1:length(t)
    clf
    xlim([0,2]);
    ylim([-1,1])
    hold on
    q1 = quiver(0,0,cos(q(i,1)),sin(q(i,1)), ShowArrowHead='off', AutoScale='off');
    q2 = quiver(cos(q(i,1)),sin(q(i,1)),cos(q(i,1))+cos(q(i,1)+q(i,2)),sin(q(i,1))+sin(q(i,1)+q(i,2)), ShowArrowHead='off', AutoScale='off');
    hold off
    pause(t(end)/length(t))
end

function dq_dt = dqdt(t,q)
    v = (15*cos(q(1) + q(2))^2 - 33*cos(q(1)) - 33*cos(q(1) + q(2)) - 2*sin(q(1))*sin(t) + 15*cos(q(1))^2 + 3*q(3)^2*cos(q(1) + q(2))^2 + 3*q(4)^2*cos(q(1) + q(2))^2 + 3*q(3)^2*sin(q(1) + q(2))^2 + 3*q(4)^2*sin(q(1) + q(2))^2 + 3*q(3)^2*cos(q(1))^2 + 3*q(3)^2*sin(q(1))^2 + 30*cos(q(1) + q(2))*cos(q(1)) - 2*sin(q(1) + q(2))*sin(t) + 3*q(3)*cos(q(1) + q(2))*sin(q(1) + q(2)) + 3*q(4)*cos(q(1) + q(2))*sin(q(1) + q(2)) + 6*q(3)*q(4)*cos(q(1) + q(2))^2 + 6*q(3)*q(4)*sin(q(1) + q(2))^2 + 3*q(3)*cos(q(1) + q(2))*sin(q(1)) + 3*q(3)*sin(q(1) + q(2))*cos(q(1)) + 3*q(4)*sin(q(1) + q(2))*cos(q(1)) + 3*q(3)*cos(q(1))*sin(q(1)) + 6*q(3)^2*cos(q(1) + q(2))*cos(q(1)) + 3*q(4)^2*cos(q(1) + q(2))*cos(q(1)) + 6*q(3)^2*sin(q(1) + q(2))*sin(q(1)) + 3*q(4)^2*sin(q(1) + q(2))*sin(q(1)) + 6*q(3)*q(4)*cos(q(1) + q(2))*cos(q(1)) + 6*q(3)*q(4)*sin(q(1) + q(2))*sin(q(1)))/(3*(cos(q(1) + q(2))*sin(q(1)) - sin(q(1) + q(2))*cos(q(1))));
    dq_dt(1:2,1) = q(3:4);
    dq_dt(3:4,1) = v;
end