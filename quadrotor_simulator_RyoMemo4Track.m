%% Initializations
clc; clear; close all;

g = 9.81;   % The gravitational acceleration [m/s^2]
l =  0.2;   % Distance from the center of mass to each rotor [m]
m =  0.5;   % Total mass of the quadrotor [kg]
I = [1.24, 1.24, 2.48]; % Mass moment of inertia [kg m^2]
mu = 3.0;   % Maximum thrust of each rotor [N]
sigma = 0.01; % The proportionality constant relating thrust to torque [m]


p = [g l m I mu sigma];

r = [0; 0; 0];
n = [0; 0; 0];


[A,B,K] = quadrotor_modeling;
u = @(z, zd) p(3) * p(1) / 4 + K*(zd - z);
% u = @(t,z,p,zd, A, B, K) quadrotor_feedback_linearization(t,z,p,zd, A,B,K);

%% Solving the initial-value problem
max = 2*pi;
speed = 100;
% tspan = linspace(0, speed*max, speed*max*100);
tspan = linspace(0, 300, 3000);
% Initial conditions
z0 = zeros(12,1);

x_ini = 0;
y_ini = 5;
z_ini = 5;
x_d_ini = 0;
y_d_ini = 0;
z_d_ini = 0;
position_init = [x_ini; y_ini; z_ini];
velocity_init = [x_d_ini; y_d_ini; z_d_ini];
t1 = 100; % time when target enters the area
t2 = 110; % time when target acceleration stops

x_d_const = 1/30;
a_vel = x_d_const/(t2-t1);
b_vel = -x_d_const*t1/(t2-t1);

x_vel = @(t) x_d_ini + (a_vel*t + b_vel)*heaviside(t-t1) - ((a_vel*t + b_vel)*heaviside(t-t1))*heaviside(t-t2) + (x_d_const)*heaviside(t-t2); 
x_pos = @(t) x_ini + (a_vel*t + b_vel)*heaviside(t-t1)*(t-t1) - (a_vel*t + b_vel)*heaviside(t-t2)*(t-t2) + (x_d_const)*heaviside(t-t2)*(t-t2);

% x_pos = @(t) x_ini;
% x_vel = @(t) x_d_ini;

%
zd = @(t) [
            % x_ini + x_d_const*heaviside(t-t1)*(t-t1);
            % x_ini + (a_vel*t + b_vel)*heaviside(t-t1) + (a_vel*t1 + b_vel)*heaviside(t-t2)*(t-t2);
            x_pos(t) - 10*heaviside(t-t1);
%             position_init;
%            y_ini + 0 * heaviside(t-t1);
            -x_pos(t);
           -x_pos(t) + 10*heaviside(t-t1);
           zeros(3,1);
%            z_ini + 0 * heaviside(t-t1);
        %    x_d_ini + x_d_const * heaviside(t-t1);
           x_vel(t);
%            velocity_init;
            -x_vel(t);
            -x_vel(t);
%            y_d_ini + 0 * heaviside(t-t1);           
%            z_d_ini + 0 * heaviside(t-t1);
           zeros(3,1)];
%}

% zd = @(t) [10*sin(t/speed); 10*cos(t/speed); sin(t/(speed/10))+10*t/tspan(end);zeros(9,1)];
% Linear form
[t,z] = ode45(@(t,z) quadrotor(t, z, u(z, zd(t)), p, r, n), tspan, z0);
% [t2,z2] = ode45(@(t,z) quadrotor(t, z, u(z, z0), p, r, n), tspan, z(end,:));
% Feedback Linearization
% [t,z] = ode45(@(t,z) quadrotor(t, z, u(t,z, p, zd(t), A,B,K), p, r, n), tspan, z0);
% [t2,z2] = ode45(@(t,z) quadrotor(t, z, u(t,z, p, z0, A,B,K), p, r, n), tspan, z(end,:));
% t = [t;t2+t(end)];
% z = [z;z2];

%% Plotting the results

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

% return
%%



%% trajectory vs measured


trajectory = zeros(length(tspan),12);
idx = 1;
    for time = tspan
        trajectory(idx,:) = zd(time);
        idx=idx+1;
    end


%

t_start = 1;
t_end = 250;
t_res = 40;
 for t_count=t_start*10:t_res:(t_end)*10
    figure(2)
    set(gcf, 'Position',  [100, 100, 1000, 400])
    subplot(1,2,1)
    hold on
    % hold on
    scatter3(z(t_count,1), z(t_count,2), z(t_count,3),'filled','MarkerFaceColor','blue');
    hold on
    scatter3(trajectory(t_count,1),trajectory(t_count,2),trajectory(t_count,3),'filled','MarkerFaceColor','red');
    view(2)
    % hold off
    xlabel("x1")
    ylabel("x2")
    zlabel("x3")
    xlim([-10 10])
    ylim([-10 10])
    zlim([-10 10])
    legend("measured", "desired")
    grid on
    disp(t_count)
    title(sprintf('t = %.1f', t(t_count)));
    % axis equal
    % drawnow
    % pause(0.01)

    subplot(1,2,2)
    hold on
    % hold on
    scatter3(z(t_count,1), z(t_count,2), z(t_count,3),'filled','MarkerFaceColor','blue');
    hold on
    scatter3(trajectory(t_count,1),trajectory(t_count,2),trajectory(t_count,3),'filled','MarkerFaceColor','red');
    view(3)
    % hold off
    xlabel("x1")
    ylabel("x2")
    zlabel("x3")
    xlim([-10 10])
    ylim([-10 10])
    zlim([-10 10])
    % legend("measured", "desired")
    grid on
    disp(t_count)
    title(sprintf('t = %.1f', t(t_count)));
    % axis equal
    drawnow
    pause(0.01)

end 
%}



figure
error = z(:,1:3) - trajectory(:,1:3);
error_vel = z(:,7:9) - trajectory(:,7:9);
subplot(2,1,1)
hold on
plot(t,error(:,1),"LineWidth",1.5);
plot(t,error(:,2),"LineWidth",1.5);
plot(t,error(:,3),"LineWidth",1.5);
legend("x1","x2","x3")
ylabel("error position [m]")
xlabel("time [s]")
grid on

subplot(2,1,2)
hold on
plot(t,error_vel(:,1),"LineWidth",1.5);
plot(t,error_vel(:,2),"LineWidth",1.5);
plot(t,error_vel(:,3),"LineWidth",1.5);
legend("x_{d1}", "x_{d2}", "x_{d3}")
ylabel("error velocity [m/s]")
% xlabel("time [s]")
grid on

% axis equal

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