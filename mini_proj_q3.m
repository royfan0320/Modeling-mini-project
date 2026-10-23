clear; clc; close all;

%% 1. Parameters
% Motor and mechanics
p.rm = 0.011;
p.Lm = 0.01;
p.Rm = 90;
p.Jm = 8e-6;
p.Km = 0.07;
p.m  = 0.06;
p.mu_r = 0.003;
p.mu_x = 0.002;
p.mu_y = 0.001;
p.bx = p.mu_x + 2*p.mu_r;
p.by = p.mu_y + 2*p.mu_r;
p.dx = 0.30;
p.dy = 0.30;

% Nozzle
p.Rn = 8;
p.Rtn = 9;
p.TeC = 300 - 273.15;
p.Ln = 0.1;
p.Cn = 0.02;
p.Ctn = 12;
p.Tlow = 50;
p.Thigh = 240;

% Events and motion
p.tColor = 15;
p.tFilament = Inf;    % Disable second event for the first full-cycle test
p.tWipe = 5;
p.tLoad = 8;
p.xTower = 0;
p.yTower = 0;
p.posTol = 1e-3;
p.velTol = 1e-3;
p.kp = 350;
p.kd = 35;
p.uMax = 24;

% Heater: steady-state feedforward plus proportional correction
p.uHeatMax = 24;
p.kHeat = 0.12;       % Test controller gain
p.tEnd = 650;        % Cooling from 240 to 50 degC takes ~250 s
p.maxEvents = 100;

% Discrete modes
p.PRINTING = 1;
p.STOPPING = 2;
p.MOVING = 3;
p.WIPING = 4;
p.COOLING = 5;
p.LOADING = 6;
p.HEATING = 7;
p.RETURNING = 8;
modeNames = {'Printing','Stopping','Move to tower','Wiping', ...
             'Cooling','Loading','Heating','Returning'};

%% 2. Initial conditions
% Start at a consistent heater electrical equilibrium for Tref=240 degC.
uEq = sqrt(p.Rn*(p.Thigh-p.TeC)/p.Rtn);
x0 = [0.10; 0; 0; 0.10; 0; 0; p.Thigh; uEq; uEq/p.Rn; 0.10; 0.10; 0];

q = p.PRINTING;
color = 1; % 1 = red, 2 = green
colorDone = false;
filamentDone = false;
t0 = 0;

T = t0;
X = x0.';
Q = q;
COLORS = color;
eventTimes = [];
eventLabels = {};
baseOpts = odeset('RelTol',1e-7, 'AbsTol',1e-9, 'MaxStep',0.25);

%% 3. Event-driven simulation
for k = 1:p.maxEvents
    if t0 >= p.tEnd
        break;
    end

    % Events may already be enabled at the beginning of a mode.
    % This is especially important for events scheduled while not printing.
    instantEvent = readyEvent(t0,x0,q,p,colorDone,filamentDone);
    if instantEvent > 0
        tEvent = t0;
        xEvent = x0;
        eventID = instantEvent;
    else
        opts = odeset(baseOpts,'Events', ...
            @(t,x) guard(t,x,q,p,colorDone,filamentDone));

        [ts,xs,te,xe,ie] = ode15s( ...
            @(t,x) StateChange(t,x,q,p), [t0 p.tEnd], x0, opts);

        T = [T; ts(2:end)];
        X = [X; xs(2:end,:)];
        Q = [Q; repmat(q,numel(ts)-1,1)];
        COLORS = [COLORS; repmat(color,numel(ts)-1,1)];

        if isempty(te)
            break;
        end
        tEvent = te(end);
        xEvent = xe(end,:).';
        eventID = ie(end);
    end

    [xNew,qNew,colorNew,eventLabel] = jump(xEvent,q,eventID,color,p);

    if q == p.PRINTING
        if eventID == 1, colorDone = true; end
        if eventID == 2, filamentDone = true; end
    end

    fprintf('%8.3f s: %s -> %s (%s)\n', ...
        tEvent,modeNames{q},modeNames{qNew},eventLabel);
    eventTimes(end+1,1) = tEvent;
    eventLabels{end+1,1} = eventLabel;

    T = [T;tEvent];
    X = [X;xNew.'];
    Q = [Q;qNew];
    COLORS = [COLORS;colorNew];

    if qNew == q && norm(xNew-xEvent) < 1e-12
        error('Transition did not change mode or state (possible Zeno loop).');
    end
    t0 = tEvent;
    x0 = xNew;
    q = qNew;
    color = colorNew;
end

if k == p.maxEvents && t0 < p.tEnd
    warning('Reached maxEvents before tEnd. Check for repeated events.');
end

%% 4. Results
figure('Color','w','Position',[100 100 1150 850]);
subplot(4,1,1);
plot(T,X(:,1),'b',T,X(:,4),'m','LineWidth',1.2);
ylabel('Position (m)'); legend('x','y','Location','best'); grid on;

subplot(4,1,2);
plot(T,X(:,7),'LineWidth',1.3); hold on;
yline(p.Tlow,'b--'); yline(p.Thigh,'r--');
ylabel('Nozzle (degC)'); grid on;

subplot(4,1,3);
stairs(T,Q,'LineWidth',1.4);
ylim([0.5 8.5]); yticks(1:8); yticklabels(modeNames);
ylabel('Mode'); grid on;

subplot(4,1,4);
stairs(T,COLORS,'LineWidth',1.4);
yticks([1 2]); yticklabels({'Red','Green'}); ylim([0.5 2.5]);
ylabel('Filament'); xlabel('Time (s)'); grid on;
sgtitle('Q3.1 Hybrid Printer: Filament Replacement');

figure('Color','w'); hold on; grid on; axis equal;
plot(X(:,1),X(:,4),'Color',[0.6 0.6 0.6]);
redIdx = Q==p.PRINTING & COLORS==1;
greenIdx = Q==p.PRINTING & COLORS==2;
plot(X(redIdx,1),X(redIdx,4),'r.','MarkerSize',7);
plot(X(greenIdx,1),X(greenIdx,4),'g.','MarkerSize',7);
plot(p.xTower,p.yTower,'kx','MarkerSize',12,'LineWidth',2);
xlim([-0.01 p.dx+0.01]); ylim([-0.01 p.dy+0.01]);
xlabel('x (m)'); ylabel('y (m)'); title('Tool-head trajectory');
legend('Full motion','Printing red','Printing green','Wipe tower', ...
       'Location','best');

%% Functions
function dx = StateChange(t,x,q,p)
    [ux,uy,un] = inputSignal(t,x,q,p);
    Me = p.m + 2*p.Jm/p.rm^2;

    Ax = [0, 1, 0;
          0, -p.bx/Me, p.Km/(p.rm*Me);
          0, -2*p.Km/(p.rm*p.Lm), -p.Rm/p.Lm];
    Ay = [0, 1, 0;
          0, -p.by/Me, p.Km/(p.rm*Me);
          0, -2*p.Km/(p.rm*p.Lm), -p.Rm/p.Lm];
    An = [-1/(p.Ctn*p.Rtn), 0, 0;
          0, -1/(p.Cn*p.Rn), 1/p.Cn;
          0, -1/p.Ln, 0];

    A = blkdiag(Ax,Ay,An,zeros(3));
    B = [0 0 0; 0 0 0; 1/p.Lm 0 0;
         0 0 0; 0 0 0; 0 1/p.Lm 0;
         0 0 0; 0 0 0; 0 0 1/p.Ln;
         zeros(3,3)];

    g = zeros(12,1);
    g(7) = x(8)^2/(p.Ctn*p.Rn) + p.TeC/(p.Ctn*p.Rtn);
    g(12) = double(q==p.WIPING || q==p.LOADING);

    dx = A*x + B*[ux;uy;un] + g;
end

function [ux,uy,un] = inputSignal(t,x,q,p)
    switch q
        case p.PRINTING
            xr = 0.15 + 0.04*sin(2*pi*t/12);
            yr = 0.15 + 0.04*cos(2*pi*t/12);
        case p.STOPPING
            xr = x(10);
            yr = x(11); % Hold saved print location
        case p.MOVING
            xr = p.xTower;
            yr = p.yTower;
        case p.RETURNING
            xr = x(10);
            yr = x(11);
        otherwise
            xr = x(1);
            yr = x(4);
    end

    if any(q==[p.PRINTING p.STOPPING p.MOVING p.RETURNING]) % Apply PD controller
        ux = saturate(p.kp*(xr-x(1))-p.kd*x(2),p.uMax);
        uy = saturate(p.kp*(yr-x(4))-p.kd*x(5),p.uMax);
    else
        ux = 0; uy = 0;
    end

    if q==p.COOLING || q==p.LOADING
        un = 0;
    else
        uFF = sqrt(p.Rn*max(p.Thigh-p.TeC,0)/p.Rtn);
        un = max(0,min(p.uHeatMax,uFF+p.kHeat*(p.Thigh-x(7))));
    end
end

function u = saturate(u,uMax)
    u = max(-uMax,min(uMax,u));
end

function [value,isterminal,direction] = guard(t,x,q,p,colorDone,filamentDone)
    switch q
        case p.PRINTING
            value = [1;1];
            if ~colorDone
                value(1) = t-p.tColor;
            end
            if ~filamentDone
                value(2) = t-p.tFilament;
            end
            direction = [1;1];
        case p.STOPPING
            value = max(abs([x(2),x(5)]))-p.velTol;
            direction = -1;
        case p.MOVING
            value = targetError(x,p.xTower,p.yTower,p);
            direction = -1;
        case p.WIPING
            value = x(12)-p.tWipe;
            direction = 1;
        case p.COOLING
            value = x(7)-p.Tlow;
            direction = -1;
        case p.LOADING
            value = x(12)-p.tLoad;
            direction = 1;
        case p.HEATING
            value = x(7)-p.Thigh;
            direction = 1;
        case p.RETURNING
            value = targetError(x,x(10),x(11),p);
            direction = -1;
        otherwise
            error('Unknown discrete mode: %d',q);
    end
    isterminal = ones(size(value));
end

function eventID = readyEvent(t,x,q,p,colorDone,filamentDone)
    eventID = 0;
    switch q
        case p.PRINTING
            if ~colorDone && t>=p.tColor
                eventID = 1;
            elseif ~filamentDone && t>=p.tFilament
                eventID = 2;
            end
        case p.STOPPING
            if max(abs([x(2),x(5)]))<=p.velTol, eventID=1; end
        case p.MOVING
            if targetError(x,p.xTower,p.yTower,p)<=0, eventID=1; end
        case p.WIPING
            if x(12)>=p.tWipe, eventID=1; end
        case p.COOLING
            if x(7)<=p.Tlow, eventID=1; end
        case p.LOADING
            if x(12)>=p.tLoad, eventID=1; end
        case p.HEATING
            if x(7)>=p.Thigh, eventID=1; end
        case p.RETURNING
            if targetError(x,x(10),x(11),p)<=0, eventID=1; end
    end
end

function e = targetError(x,xr,yr,p)
    e = max([abs(x(1)-xr)/p.posTol, abs(x(4)-yr)/p.posTol, ...
             abs(x(2))/p.velTol, abs(x(5))/p.velTol])-1;
end

function [xn,qn,cn,label] = jump(x,q,event,color,p)
    xn = x; cn = color;
    switch q
        case p.PRINTING
            qn = p.STOPPING;
            xn(10) = x(1); xn(11) = x(4);
            if event==1
                label = 'colour change requested';
            else
                label = 'filament finished';
            end
        case p.STOPPING
            qn = p.MOVING; label = 'head stopped';
        case p.MOVING
            qn = p.WIPING; xn(12) = 0; label = 'arrived at wipe tower';
        case p.WIPING
            qn = p.COOLING; label = 'wiping completed';
        case p.COOLING
            qn = p.LOADING; xn(12) = 0; label = 'nozzle cooled';
        case p.LOADING
            qn = p.HEATING; cn = 2; label = 'new filament loaded';
        case p.HEATING
            qn = p.RETURNING; label = 'nozzle reheated';
        case p.RETURNING
            qn = p.PRINTING; label = 'returned to saved position';
        otherwise
            error('Unknown discrete mode: %d',q);
    end
end
