clear; clc; close all;

% Parameters
p.rm=0.010;
p.Lm=0.004;
p.Rm=2.0;
p.Jm=1.2e-5;
p.Km=0.050;
p.m=0.120;
p.dx=0.30;
p.dy=0.30;
p.mu_r=0.010;
p.mu_x=0.020;
p.mu_y=0.020;
p.bx=8*p.mu_r+2*p.mu_x; p.by=4*p.mu_r+2*p.mu_y;
p.Rn=8;
p.Rtn=9;
p.Te=25;
p.Ln=0.1;
p.Cn=0.02;
p.Ctn=6;
p.Tlow=50;
p.Thigh=240;
p.tempTol=0.5;
p.tColor=4.8;
p.tFilament=203.5;
p.tWipe=5;
p.tLoad=8;

p.posTol=1e-3;
p.velTol=1e-3;
p.release=1e-4;
p.xTower=0;
p.yTower=0;
p.homeInput=[-5;-5];
p.sampleTime=0.02;
p.startupTime=0.10;
p.kp=350;
p.kd=35;
p.kpSquare=60;
p.kdSquare=15;
p.uMax=24;
p.kHeat=0.12;
p.corners=[.10 .10; .20 .10; .20 .20; .10 .20; .10 .10];
p.tEnd=450;
p.maxEvents=20000;

% Initialize and simulate
uEq=sqrt(p.Rn*(p.Thigh-p.Te)/p.Rtn);
x0=[p.corners(1,1);0;0;p.corners(1,2);0;0;p.Thigh;uEq;uEq/p.Rn;.10;.10;0];

d=struct('q',1,'v',0,'s',1,'color',1,'done',[false false],'wall',false(4,1),'nextSample',Inf,'homeStart',0);

[T,X,D,events,d]=simulate(x0,d,p);

% Results
modeNames={'Printing','Stopping','Move to tower','Wiping','Cooling','Loading','Heating','Returning'};

figure('Color','w');
subplot(4,1,1);
plot(T,X(:,1),T,X(:,4));
ylabel('Position (m)');
legend('x','y');
grid on;

subplot(4,1,2);
plot(T,X(:,7));
yline(p.Tlow,'--');
yline(p.Thigh,'--');
ylabel('Nozzle (degC)');
grid on;

subplot(4,1,3);
stairs(T,D(:,1));
yticks(1:8);
yticklabels(modeNames);
ylabel('Main mode');
grid on;

subplot(4,1,4);
stairs(T,D(:,4));
yticks([1 2]);
yticklabels({'Red','Green'});
ylabel('Color');
xlabel('Time (s)');
grid on;


figure('Color','w');
hold on;
axis equal;
grid on;

plot(p.corners(:,1),p.corners(:,2),'k--');
plot(X(:,1),X(:,4),'Color',[.65 .65 .65]);
r=D(:,1)==1 & D(:,4)==1;
g=D(:,1)==1 & D(:,4)==2;
plot(X(r,1),X(r,4),'r.');
plot(X(g,1),X(g,4),'g.');
plot(0,0,'kx');

xlim([-.01 p.dx+.01]);
ylim([-.01 p.dy+.01]);
xlabel('x (m)');
ylabel('y (m)');
legend('Desired square','Full motion','Printing red','Printing green','Wipe tower','Location','best');

figure('Color','w');
subplot(2,1,1);
stairs(T,D(:,2));
yticks(0:3);
ylabel('Q2 substate');
grid on;

subplot(2,1,2);
stairs(T,D(:,3));
yticks(1:5);
ylabel('Square state');
xlabel('Time (s)');
grid on;

fprintf('Final time %.3f s; square state %d; replacement requests processed %d/2.\n',T(end),D(end,3),sum(d.done));

% Final discrete state d is returned for inspection.
if D(end,3)~=5
    warning('The square has not finished within the simulation horizon.');
end

%% Functions
% Continuous physical model
function dx=StateChange(~,x,d,p)

    u=inputSignal(x,d,p);

    Me=p.m+2*p.Jm/p.rm^2;
    Ax=[0 1 0; 0 -p.bx/Me p.Km/(p.rm*Me); 0 -2*p.Km/(p.rm*p.Lm) -p.Rm/p.Lm];
    Ay=Ax;
    Ay(2,2)=-p.by/Me;
    An=[-1/(p.Ctn*p.Rtn) 0 0; 0 -1/(p.Cn*p.Rn) 1/p.Cn; 0 -1/p.Ln 0];

    A=blkdiag(Ax,Ay,An,zeros(3));
    
    B=zeros(12,3);
    B(3,1)=1/p.Lm;
    B(6,2)=1/p.Lm;
    B(9,3)=1/p.Ln;

    g=zeros(12,1);
    g(7)=x(8)^2/(p.Ctn*p.Rn)+p.Te/(p.Ctn*p.Rtn);
    g(12)=double(d.q==4 || d.q==6);

    free=A*x+B*u+g;

    blockX=(d.wall(1) && x(1)<=1e-10 && x(2)<=0 && free(2)<0) || (d.wall(2) && x(1)>=p.dx-1e-10 && x(2)>=0 && free(2)>0);
    blockY=(d.wall(3) && x(4)<=1e-10 && x(5)<=0 && free(5)<0) || (d.wall(4) && x(4)>=p.dy-1e-10 && x(5)>=0 && free(5)>0);

    rows=[];
    if blockX
        rows=[rows 1 2];
    end

    if blockY
        rows=[rows 4 5];
    end

    A(rows,:)=0;
    B(rows,:)=0;
    g(rows)=0;
    dx=A*x+B*u+g;
end

% Inputs
function u=inputSignal(x,d,p)

    ux=0;
    uy=0;
    
    switch d.q
        case 1
            target=p.corners(min(d.s+1,5),:).';
            kp=p.kpSquare; kd=p.kdSquare;
            uxy=kp*(target-x([1 4]))-kd*x([2 5]);
            ux=uxy(1);
            uy=uxy(2);

        case {2,8}
            uxy=p.kp*(x([10 11])-x([1 4]))-p.kd*x([2 5]);
            ux=uxy(1);
            uy=uxy(2);

        case 3
            ux=p.homeInput(1);
            uy=p.homeInput(2);
            if d.v==1 || d.v==3
                ux=0;
            end
            if d.v==2 || d.v==3
                uy=0;
            end
    end

    un=0;

    if d.q~=5 && d.q~=6
        uFF=sqrt(p.Rn*(p.Thigh-p.Te)/p.Rtn);
        un=uFF+p.kHeat*(p.Thigh-x(7));
    end
    u=[max(-p.uMax,min(p.uMax,[ux;uy]));max(0,min(p.uMax,un))];
end

% Guards
function [value,terminal,direction]=guard(t,x,d,p)

    walls=[x(1);x(1)-p.dx;x(4);x(4)-p.dy];
    walls(d.wall)=1;

    release=[x(1);p.dx-x(1);x(4);p.dy-x(4)]-p.release;
    release(~d.wall)=1;
    value=[walls;release;ones(4,1)];

    direction=[-1;1;-1;1;ones(6,1);-1;-1];
    terminal=ones(12,1);

    if d.q==1
        if ~d.done(1) && isfinite(p.tColor)
            value(9)=t-p.tColor;
        end
        if ~d.done(2) && isfinite(p.tFilament)
            value(10)=t-p.tFilament;
        end
        if d.s<=4
            err=max(abs(x([1 4])-p.corners(d.s+1,:).'))/p.posTol;
            speed=max(abs(x([2 5])))/p.velTol;
            value(12)=max(err,speed)-1;
        end
    end

    switch d.q
        case 2
            value(11)=max(abs(x([2 5])))-p.velTol;

        case 3
            value(11)=t-d.nextSample;
            direction(11)=1;

        case 4
            value(11)=x(12)-p.tWipe;
            direction(11)=1;

        case 5
            value(11)=x(7)-p.Tlow;

        case 6
            value(11)=x(12)-p.tLoad;
            direction(11)=1;

        case 7
            value(11)=x(7)-(p.Thigh-p.tempTol);
            direction(11)=1;

        case 8
            err=max(abs(x([1 4])-x([10 11])))/p.posTol;
            speed=max(abs(x([2 5])))/p.velTol;
            value(11)=max(err,speed)-1;
    end
end

% Jumps
function [x,d,label]=jump(t,x,d,id,p)

    label='';

    if id<=4
        if id==1
            x(1)=0; x(2)=0;
        end

        if id==2
            x(1)=p.dx; x(2)=0;
        end

        if id==3
            x(4)=0; x(5)=0;
        end

        if id==4
            x(4)=p.dy; x(5)=0;
        end

        d.wall(id)=true; label=sprintf('wall %d: velocity reset; main mode unchanged',id);
        return

    elseif id<=8
        d.wall(id-4)=false; % allow another collision after moving inward
        return

    elseif id==9 || id==10
        d.done(id-8)=true;
        d.q=2;
        x([10 11])=x([1 4]);
        label=sprintf('request %d: save position and stop',id-8);
        return

    elseif id==12
        d.s=min(d.s+1,5);
        label=sprintf('square state -> %d',d.s);
        return
    end

    switch d.q
        case 2
            d.q=3;
            d.v=0;
            d.homeStart=t;
            d.nextSample=t+p.sampleTime;
            label='start fixed-input homing';

        case 3
            % Motor sensor model: omegaL=(vx+vy)/rm, omegaR=(vx-vy)/rm
            omega=[x(2)+x(5);x(2)-x(5)]/p.rm;
            velocity=(p.rm/2)*[1 1;1 -1]*omega;
            old=d.v;

            if t-d.homeStart>=p.startupTime-1e-10
                % + and 0 mean moving and stopped along the commanded axis
                if d.v~=1 && d.v~=3 && abs(velocity(1))<=p.velTol
                    d.v=bitor(d.v,1); % q0->q1 or q2->q3 on (H,0)
                end
                if d.v~=2 && d.v~=3 && abs(velocity(2))<=p.velTol
                    d.v=bitor(d.v,2); % q0->q2 or q1->q3 on (V,0)
                end
            end

            d.nextSample=d.nextSample+p.sampleTime;

            if d.v~=old
                label=sprintf('Q2 state q%d -> q%d from measured speed',old,d.v);
            end

            if d.v==3
                d.q=4;
                x(12)=0;
                d.nextSample=Inf;
                label='Q2 q3 reached: homing complete; start wiping';
            end
            
        case 4
            d.q=5;
            label='wiping finished; heater off';

        case 5
            d.q=6;
            x(12)=0;
            label='cooled; start loading';

        case 6
            d.q=7;
            d.color=2;
            label='loaded; reheat';

        case 7
            d.q=8;
            label='heated; return to saved position';

        case 8
            d.q=1;
            label='saved position AND low speed reached; resume printing';

    end
end

% Simulate
function [T,X,D,events,d]=simulate(x,d,p)

    t=0;
    T=t;
    X=x.';
    D=[d.q d.v d.s d.color];
    events=cell(0,3);
    options=odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',0.1);

    for iteration=1:p.maxEvents
        if t>=p.tEnd
            break;
        end

        [v,~,directions]=guard(t,x,d,p);
        active=false(12,1);
        active(9:10)=d.q==1 & ~d.done(:) & isfinite([p.tColor;p.tFilament]) & v(9:10)>=-1e-10;
        active(12)=d.q==1 && d.s<=4 && v(12)<=1e-10;
        active(11)=v(11)<=1e-10 && d.q~=1;

        if directions(11)==1
            active(11)=v(11)>=-1e-10;
        end

        u=inputSignal(x,d,p);
        active(1)=~d.wall(1) && x(1)<=1e-10 && (x(2)<0 || (x(2)<=0 && u(1)<0));
        active(2)=~d.wall(2) && x(1)>=p.dx-1e-10 && (x(2)>0 || (x(2)>=0 && u(1)>0));
        active(3)=~d.wall(3) && x(4)<=1e-10 && (x(5)<0 || (x(5)<=0 && u(2)<0));
        active(4)=~d.wall(4) && x(4)>=p.dy-1e-10 && (x(5)>0 || (x(5)>=0 && u(2)>0));
        hits=find(active);

        if isempty(hits)
            opts=odeset(options,'Events',@(tt,xx)guard(tt,xx,d,p));
            [ts,xs,te,xe,ie]=ode15s(@(tt,xx)StateChange(tt,xx,d,p),[t p.tEnd],x,opts);
            
            if isempty(te)
                T=[T;ts(2:end)];
                X=[X;xs(2:end,:)];
                D=[D;repmat([d.q d.v d.s d.color],numel(ts)-1,1)];
                t=ts(end);
                break
            end

            next=min(te);
            same=find(abs(te-next)<=1e-9);
            hits=sort(unique(ie(same)));
            x=xe(same(1),:).';
            keep=find(ts>t & ts<next);
            T=[T;ts(keep);next];
            X=[X;xs(keep,:);x.'];
            D=[D;repmat([d.q d.v d.s d.color],numel(keep)+1,1)];
            t=next;
        end

        walls=hits(hits<=4);

        if isempty(walls)
            selected=hits(1);
        else
            selected=walls;
        end

        for id=selected(:).'

            before=d.q;
            [x,d,label]=jump(t,x,d,id,p);

            if ~isempty(label)
                events(end+1,:)={t,id,label};
                fprintf('%9.3f s: mode %d -> %d | %s\n',t,before,d.q,label);
            end
        end

        T=[T;t];
        X=[X;x.'];
        D=[D;d.q d.v d.s d.color];

    end

    if t<p.tEnd-1e-8
        error('Event limit reached before tEnd; check repeated events.');
    end
end
