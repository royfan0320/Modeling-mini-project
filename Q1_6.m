%% Q1.6
clear; clc; close all;

%% Assumed parameters
rho = 1250; % kg/m^3
c_p = 1800; % J/(kg*K)
k = 0.20; % W/(m*K)
width = 0.0004; % Square cross-section side length, m
h_air = 10; % W/(m^2*K)
h_bed = 100; % W/(m^2*K)
T_e = 25; % Air temperature, deg C
T_b = 60; % Bed temperature, deg C
T_dep = 210; % Deposition temperature, deg C
tol = 1; % Allowed temperature difference around the layer, K

alpha = k/(rho*c_p);
beta = (h_bed+3*h_air)/(rho*c_p*width);
T_eq = (h_bed*T_b+3*h_air*T_e)/(h_bed+3*h_air);

%% Ellipse and length
a = 0.08; b = 0.05;
t_print = 12;
w = 2*pi/t_print;
tg = linspace(0,t_print,10001);
speed = sqrt((a*w*sin(w*tg)).^2 + (b*w*cos(w*tg)).^2);
sg = cumtrapz(tg,speed);
L = sg(end);

N = 200;
ds = L/N;
s = (0:N-1)'*ds;
dt = 0.01;
t = 0:dt:80;
printed_length = interp1(tg,sg,min(t,t_print));
if dt*(2*alpha/ds^2+beta) > 1
    error('Use a smaller time step.');
end

%% Temperature calculation
T = nan(N,1); % Unprinted material has no temperature
T_history = nan(N,length(t));
T_max = zeros(size(t));
T_min = zeros(size(t));
old_active = 0;
t_uniform = NaN;

for it = 1:length(t)
    active = sum(s <= printed_length(it)+1e-12); % Allow for round-off
    T(old_active+1:active) = T_dep;
    old_active = active;
    if t(it) < t_print
        T(active) = T_dep;
    end

    T_history(:,it) = T;
    T_max(it) = max(T(1:active));
    T_min(it) = min(T(1:active));
    if t(it) >= t_print && T_max(it)-T_min(it) <= tol
        t_uniform = t(it);
        break;
    end

    T_new = T;
    if t(it) < t_print
        for j = 1:active-1
            if j == 1
                d2T = 2*(T(2)-T(1))/ds^2;
            else
                d2T = (T(j+1)-2*T(j)+T(j-1))/ds^2;
            end
            T_new(j) = T(j)+dt*(alpha*d2T-beta*(T(j)-T_eq));
        end
    else
        for j = 1:N
            left = j-1;
            right = j+1;
            if j == 1, left = N; end
            if j == N, right = 1; end
            d2T = (T(right)-2*T(j)+T(left))/ds^2;
            T_new(j) = T(j)+dt*(alpha*d2T-beta*(T(j)-T_eq));
        end
    end
    T = T_new;
end

if isnan(t_uniform)
    error('Increase the final simulation time.');
end
t = t(1:it);
T_history = T_history(:,1:it);
T_max = T_max(1:it);
T_min = T_min(1:it);
fprintf('Nearly uniform at %.2f s from the start.\n',t_uniform);
fprintf('Waiting time after printing: %.2f s.\n',t_uniform-t_print);

%% Plots

subplot(1,2,1);
img = imagesc(t,s,T_history);
set(img,'AlphaData',~isnan(T_history));
set(gca,'YDir','normal','Color',[0.9 0.9 0.9]);
colorbar; xlabel('Time (s)'); ylabel('Arc length (m)');
title('Temperature (deg C)');

subplot(1,2,2);
plot(t,T_max-T_min); grid on;
yline(tol,'--');
xlabel('Time (s)'); ylabel('Temperature difference (K)');
