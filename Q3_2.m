%% Q3.2
clear; clc; close all;

%% Same assumed parameters as Q1.4
r_m = 0.010; L_m = 0.004; R_m = 2.0;
J_m = 1.2e-5; K_m = 0.050; m = 0.120;
mu_r = 0.010; mu_l_x = 0.020; mu_l_y = 0.020;
b_x = 8*mu_r+2*mu_l_x;
b_y = 4*mu_r+2*mu_l_y;
I_theta = J_m+m*r_m^2/2;
d = r_m^2*(b_x+b_y)/4;
c = r_m^2*(b_x-b_y)/4;
A = [0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 -d/I_theta -c/I_theta K_m/I_theta 0;
     0 0 -c/I_theta -d/I_theta 0 K_m/I_theta;
     0 0 -K_m/L_m 0 -R_m/L_m 0;
     0 0 0 -K_m/L_m 0 -R_m/L_m];
B = [0 0;0 0;0 0;0 0;1/L_m 0;0 1/L_m];
H = r_m/2*[1 1;1 -1];

%% Square and controller settings
x0 = 0.08; y0 = 0.08; side = 0.10;
corners = [x0 y0; x0+side y0; x0+side y0+side; ...
           x0 y0+side; x0 y0];
Kp = 60;                 % voltage position feedback gain, V/m
Kd = 2;                  % voltage velocity feedback gain, V*s/m
position_tol = 1e-7;      % m
velocity_tol = 1e-6;      % m/s
dt = 0.00005;
t = (0:dt:30)';
X = [(x0+y0)/r_m;(x0-y0)/r_m;0;0;0;0];
position = [x0;y0];
mode = 1; % 1 right, 2 up, 3 left, 4 down, 5 done
finish_time = NaN;
position_history = zeros(length(t),2);
mode_history = zeros(length(t),1);

%% hybrid controller and motor model
for it = 1:length(t)
    omega = X(3:4); % Ideal motor angular velocity sensor outputs
    velocity = H*omega;

    if mode <= 4
        target = corners(mode+1,:)';
        if max(abs(target-position)) <= position_tol && max(abs(velocity)) <= velocity_tol
            mode = mode+1;

            if mode == 5, finish_time = t(it); end
        end
    end

    % Hold the next corner as the target until position and speed are small.
    if mode <= 4
        ref = corners(mode+1,:)';
    else
        ref = corners(5,:)';
    end
    u_axis = Kp*(ref-position)-Kd*velocity;
    u_motor = [(u_axis(1)+u_axis(2))/2; (u_axis(1)-u_axis(2))/2];

    position_history(it,:) = (H*X(1:2))';
    mode_history(it) = mode;
    if mode == 5 && t(it) >= finish_time+0.3
        break;
    end

    position = position+dt*velocity;
    X = X+dt*(A*X+B*u_motor);   % Forward Euler
end

if mode ~= 5, error('Square is not completed'); end
t = t(1:it);
position_history = position_history(1:it,:);
mode_history = mode_history(1:it);

%% Plots
subplot(1,2,1);
plot(corners(:,1),corners(:,2),'r--', ...
    position_history(:,1),position_history(:,2),'b-');
axis equal; grid on;
xlim([x0-0.01 x0+side+0.01]); ylim([y0-0.01 y0+side+0.01]);
xlabel('x (m)'); ylabel('y (m)'); title('Square path');
legend('Ideal square','Simulation','Location','best');

subplot(1,2,2);
stairs(t,mode_history); grid on;
yticks(1:5); yticklabels({'Right','Up','Left','Down','Done'});
xlabel('Time (s)'); ylabel('Discrete state');
