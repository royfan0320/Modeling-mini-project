%% Q1.4 - Ellipse and Lissajous curve
clear; clc; close all;

%% assumed parameters
r_m = 0.010;% m
L_m = 0.004;% H
R_m = 2.0;% ohm
J_m = 1.2e-5;% kg*m^2
K_m = 0.050;% N*m/A, also V*s/rad
m = 0.120;% kg
mu_r = 0.010;% N*s/m, per rolling bearing
mu_l_x = 0.020;% N*s/m, per x linear bearing
mu_l_y = 0.020;% N*s/m, per y linear bearing

b_x = 8*mu_r + 2*mu_l_x;
b_y = 4*mu_r + 2*mu_l_y;
I_theta = J_m + m*r_m^2/2;
M_e = m + 2*J_m/r_m^2;
d = r_m^2*(b_x+b_y)/4;
c = r_m^2*(b_x-b_y)/4;

%% from Q1.2
A = [0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 -d/I_theta -c/I_theta K_m/I_theta 0;
     0 0 -c/I_theta -d/I_theta 0 K_m/I_theta;
     0 0 -K_m/L_m 0 -R_m/L_m 0;
     0 0 0 -K_m/L_m 0 -R_m/L_m];
B = [0 0; 0 0; 0 0; 0 0; 1/L_m 0; 0 1/L_m];

%% Draw two shapes
max_error = zeros(1,2);
for shape = 1:2
    if shape == 1
        shape_name = 'ellipse';
        T = 12;
        a = 0.08;
        b = 0.05;
        freq_x = 1;
        freq_y = 1;
        phase_x = pi/2;
    else
        shape_name = 'Lissajous';
        T = 20;
        a = 0.060;
        b = a;
        freq_x = 2;
        freq_y = 1;
        phase_fraction = 0.369;
        phase_x = 2*pi*phase_fraction;
    end
    t = (0:0.01:T)';
    w = 2*pi/T;
    phi = w*t;
    w_x = freq_x*w;
    w_y = freq_y*w;
    x_c = 0.15;
    y_c = 0.15;

    %% Reference position and first three derivatives
    x_ref = x_c + a*sin(freq_x*phi + phase_x);
    y_ref = y_c + b*sin(freq_y*phi);
    dx = a*w_x*cos(freq_x*phi + phase_x);
    dy = b*w_y*cos(freq_y*phi);
    ddx = -a*w_x^2*sin(freq_x*phi + phase_x);
    ddy = -b*w_y^2*sin(freq_y*phi);
    dddx = -a*w_x^3*cos(freq_x*phi + phase_x);
    dddy = -b*w_y^3*cos(freq_y*phi);

    %% Input design from Q1.3
    i_x = r_m/K_m*(M_e*ddx + b_x*dx);
    i_y = r_m/K_m*(M_e*ddy + b_y*dy);
    di_x = r_m/K_m*(M_e*dddx + b_x*ddx);
    di_y = r_m/K_m*(M_e*dddy + b_y*ddy);
    u_x = L_m*di_x + R_m*i_x + 2*K_m/r_m*dx;
    u_y = L_m*di_y + R_m*i_y + 2*K_m/r_m*dy;
    u_m_l = (u_x+u_y)/2;
    u_m_r = (u_x-u_y)/2;
    u_motor = [u_m_l u_m_r];

    %% Initial state
    X0 = [(x_ref(1)+y_ref(1))/r_m;
          (x_ref(1)-y_ref(1))/r_m;
          (dx(1)+dy(1))/r_m;
          (dx(1)-dy(1))/r_m;
          (i_x(1)+i_y(1))/2;
          (i_x(1)-i_y(1))/2];

    %% Simulate with the designed voltages
    options = odeset('RelTol',1e-7,'AbsTol',1e-9);
    [~,X] = ode15s(@(time,state) A*state + ...
        B*interp1(t,u_motor,time,'linear','extrap')',t,X0,options);
    x_sim = r_m/2*(X(:,1)+X(:,2));
    y_sim = r_m/2*(X(:,1)-X(:,2));
    error_x = x_sim-x_ref;
    error_y = y_sim-y_ref;
    max_error(shape) = max(sqrt(error_x.^2+error_y.^2));
    fprintf('%s: maximum tracking error = %.3e m\n',shape_name,max_error(shape));

    %% plots
    figure;
    if isprop(gcf,'Theme'), set(gcf,'Theme','light'); end
    subplot(2,2,1);
    plot(x_ref,y_ref,'r--',x_sim,y_sim,'b-');
    axis equal; grid on;
    xlabel('x (m)'); ylabel('y (m)'); title(shape_name);
    legend('Reference','Simulation','Location','best');

    subplot(2,2,2);
    plot(t,u_m_l,t,u_m_r); grid on;
    xlabel('Time (s)'); ylabel('Voltage (V)');
    legend('u_{m,l}','u_{m,r}','Location','best');

    subplot(2,2,3);
    plot(t,x_sim,t,y_sim); grid on;
    xlabel('Time (s)'); ylabel('Position (m)');
    legend('x','y','Location','best');

    subplot(2,2,4);
    plot(t,1e6*error_x,t,1e6*error_y); grid on;
    xlabel('Time (s)'); ylabel('Error (micrometres)');
    legend('e_x','e_y','Location','best');
end
