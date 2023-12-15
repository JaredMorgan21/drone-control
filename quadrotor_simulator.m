%% Initializations
clc; clear; close all;

g = 9.81;   % The gravitational acceleration [m/s^2]
l =  0.2;   % Distance from the center of mass to each rotor [m]
m =  0.5;   % Total mass of the quadrotor [kg]
I = [1.24, 1.24, 2.48]; % Mass moment of inertia [kg m^2]
mu = 3.0;   % Maximum thrust of each rotor [N]
sigma = 0.01; % The proportionality constant relating thrust to torque [m]


p = [g l m I mu sigma];

r1 = [0; 0; 0];
n1 = [0; 0; 0];

r2 = [0; 0; 0];
n2 = [0; 0; 0];

[A,B,K] = quadrotor_modeling();
Ki = K;
Ki(:,1:3) = 0.1*Ki(:,1:3);
Ki(:,4:6) = 0*Ki(:,4:6);
Ki(:,7:9) = 0.1*Ki(:,7:9);
Ki(:,10:12) = 0.1*Ki(:,10:12);

u1 = @(z, zd) p(3) * p(1) / 4 + K*(zd(1:12) - z(1:12)) - Ki(:,1:3)*z(13:15);
u2 = @(z, zd) p(3) * p(1) / 4 + K*(zd(1:12) - z(1:12)) - Ki*z(13:24);

% u = @(t,z,p,zd, A, B, K) quadrotor_feedback_linearization(t,z,p,zd, A,B,K);

%% Solving the initial-value problem
maximum = 20;
speed = 6;
tspan = linspace(0, speed*maximum, speed*maximum*100);
% Initial conditions
z0 = zeros(24,1);
UAV = @(t) [-10 + t/speed; 10*cos(t/speed); 10*sin(t/speed);zeros(21,1)];
UAVspeed = (norm(UAV(2)-UAV(1))/(tspan(2) - tspan(1)))
% zd = @(t) [10*sin(t/speed); 10*cos(t/speed); sin(t/(speed/10))+10*t/tspan(end);zeros(21,1)];
zd = UAV;
% Linear form
options = odeset('Event', @(t,z) catchDrone(t, z, zd(t), 0.1), 'RelTol', 0.1);
[t,z, te, ze] = ode45(@(t,z) quadrotor(t, z, u1(z, zd(t)), p, r1, n1, zd(t)), tspan, z0, options);
if(isempty(te))
    disp("Failed to capture")
else
    disp("Captured! Returning to home")
    options = odeset('Event', @(t,z) catchDrone(t, z, z0, 0.1), 'RelTol', 0.1);
    [t2,z2] = ode45(@(t,z, eta) quadrotor(t, z, u2(z, z0), p, r_gen(2)', n_gen(2)', z0), tspan, z(end,:), options);
end

%% trajectory vs measured
h_trajectory = figure(1);
hold on
plot3(z(:,1), z(:,2), z(:,3),'Color','red','LineWidth',1)
plot3(z2(:,1), z2(:,2), z2(:,3),'Color','blue','LineWidth',1)
trajectory = zeros(length(tspan),24);
idx = 1;
for time = tspan
    trajectory(idx,:) = zd(time);
    idx=idx+1;
end
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3), 'Color','black','LineWidth',1,'LineStyle','--');
view(3)
hold off
xlabel("x1")
ylabel("x2")
zlabel("x3")
legend('Controlled Drone','UAV', 'Location','best')
grid on
xlim([-10 10])
ylim([-10 10])
zlim([-10 10])

%% Plotting the results
figure(2)
t = [t; t2 + t(end)];
z = [z; z2];

for i=1:4
    ax(i) = subplot(2,2,i,'NextPlot','Add','Box','on','XGrid','on','YGrid','on',...
                'Xlim',[t(1), t(end)],...
                'TickLabelInterpreter','LaTeX','FontSize',14);
    xlabel('t','Interpreter','LaTeX','FontSize',14);        
end


plot(ax(1), t,z(:,1:3), 'LineWidth', 1.5);
legend(ax(1), {'$x_1$', '$x_2$', '$x_3$'},... 
    'Interpreter', 'LaTeX', 'FontSize', 14);
title(ax(1), '${\bf x}$','Interpreter','LaTeX','FontSize',14);
xlabel(ax(1), 't','Interpreter','LaTeX','FontSize',14);

plot(ax(3), t, z(:,4:6), 'LineWidth', 1.5);
legend(ax(3), {'$\phi$', '$\theta$', '$\psi$'},...
    'Interpreter', 'LaTeX', 'FontSize', 14);
title(ax(3), '\boldmath$\alpha$','Interpreter','LaTeX','FontSize',14);

plot(ax(2), t, z(:,7:9), 'LineWidth', 1.5);
legend(ax(2), {'$\dot{x}_1$', '$\dot{x}_2$', '$\dot{x}_3$'},...
    'Interpreter', 'LaTeX', 'FontSize', 14);
title(ax(2), '$\dot{\bf x}$','Interpreter','LaTeX','FontSize',14);

plot(ax(4), t, z(:,10:12), 'LineWidth', 1.5);
legend(ax(4), {'$\omega_1$', '$\omega_2$', '$\omega_3$'},...
    'Interpreter', 'LaTeX', 'FontSize', 14);
title(ax(4), '\boldmath$\omega$','Interpreter','LaTeX','FontSize',14);

% plot3(0,0,0,'filled')

% exportgraphics(h_trajectory,'fig_trajectory.pdf','ContentType','vector')

return

%% Animation
animation_fig = figure;

airspace_box_length = 4;

animation_axes = axes('Parent', animation_fig,...
    'NextPlot','add','DataAspectRatio',[1 1 1],...
    'Xlim',airspace_box_length*[-2.5 2.5],...
    'Ylim',airspace_box_length*[-2.5 2.5],...
    'Zlim',airspace_box_length*[0 2.5],...
    'box','on','Xgrid','on','Ygrid','on','Zgrid','on',...
    'TickLabelInterpreter','LaTeX','FontSize',14);

view(animation_axes, 3);

N = 10;
Q = linspace(0,2*pi,N)';
circle = 0.3*l*[cos(Q) sin(Q) zeros(N,1)];
loc = l*[1 0 0; 0 1 0; -1 0 0; 0 -1 0];


silhouette = plot3(0,0,0, '--', 'Color', 0.5*[1 1 1], 'LineWidth', 1 ,...
    'Parent', animation_axes);
body = plot3(0,0,0, 'Color',lines(1), 'LineWidth', 2,...
        'Parent', animation_axes);
for i=1:4
    rotor(i) = plot3(0,0,0, 'Color', lines(1), 'LineWidth', 2,...
        'Parent', animation_axes);
end

tic;
for k=1:length(t)
    
    R = [ cos(z(k,5))*cos(z(k,6)), sin(z(k,4))*sin(z(k,5))*cos(z(k,6)) - cos(z(k,4))*sin(z(k,6)), sin(z(k,4))*sin(z(k,6)) + cos(z(k,4))*sin(z(k,5))*cos(z(k,6));
          cos(z(k,5))*sin(z(k,6)), cos(z(k,4))*cos(z(k,6)) + sin(z(k,4))*sin(z(k,5))*sin(z(k,6)), cos(z(k,4))*sin(z(k,5))*sin(z(k,6)) - sin(z(k,4))*cos(z(k,6));
                     -sin(z(k,5)),                                 sin(z(k,4))*cos(z(k,5)),                                 cos(z(k,4))*cos(z(k,5))];
    for i=1:4
        ctr(i,:) = z(k,1:3) + loc(i,:)*R';
        pose = ones(N,1)*z(k,1:3) + (ones(N,1)*loc(i,:) + circle)*R';
        set(rotor(i), 'XData', pose(:,1), 'YData', pose(:,2),  'ZData', pose(:,3) );
         
    end
    % plot3(zd(1,:), zd(2,:), zd(3,:), Color='r');
    set(silhouette,'XData', [0, z(k,1), z(k,1), z(k,1)],...
        'YData', [0, 0, z(k,2), z(k,2)],...
        'ZData', [0, 0, 0, z(k,3)]);
    set(body, 'XData', [ctr([1 3],1); NaN; ctr([2 4],1)], ...
        'YData', [ctr([1 3],2); NaN; ctr([2 4],2)],...
        'ZData', [ctr([1 3],3); NaN; ctr([2 4],3)] );
    pause(t(k)-toc);
    pause(0.01);
end

function [value,isterminal,direction] = catchDrone(t, z, zd, epsilon)
    value = norm( zd(1:3) - z(1:3)) > epsilon;
    isterminal = 1; % = 1 -> the integration is to terminate.
    direction = 0;
end