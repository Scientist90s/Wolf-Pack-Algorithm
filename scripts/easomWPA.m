clc;
clear all;
close all;

% Easom function
easom = @(x1,x2) (-cos(x1).*cos(x2).*(exp(-((x1-pi).^2+(x2-pi).^2))));

plotfn(easom, [-8,8], [-8,8]);
plotfn(easom, [-8,8], [-8,8]);
view(0,90)

%%% Wolf pack algorithm parameters initialization
% number of wolfs
N = 20;
% number of dimensions
D = 2;
% maximum number of iterations
kmax = 30; k=1;
% step coefficient
S = 0.12;
% Distance determinant coefficient
Lnear = 0.3; %0.08
% max number of  iterations in scouting behaviour
Tmax = 8; T=0;
% population renewing proportional constant
beta = 2;
% steps
stepa = S; stepb = 2*S; stepc = S/2;
% h limits
hmin = 10; hmax= 20;
% R to delete wolves from end
Rmin = N/(2*beta); Rmax = N/beta;

%%% initial position
xmin = -5;
xmax = 5;
X = (xmax-xmin).*rand(N,D) + xmin;
m = X(:,1); n = X(:,2);
hold on
p = plot(m, n, "or");
hold off

p.XDataSource = "m";
p.YDataSource = "n";

%%% starting iterations
while k<kmax
    step = 1; overfl=0;
    while ~overfl
        if step==1
            % finding lead
            fitness = easom(X(:,1),X(:,2));
            [lead, idx_lead] = min(fitness);

            % data for plots
            [best(k), idx_best] = min(fitness);
            worst(k) = max(fitness);
            avg(k) = mean(fitness);

            step=2;
        elseif step==2
            % scouting
            brkfl = 1;
            while T<Tmax && brkfl
                for i=1:N
                    if i~=idx_lead
                        h = randi([hmin,hmax], 1);
                        movefl = 0; x1new = -6; x2new = -6; fitnessnew = fitness(i);
                        for j=1:h
                            X1new = X(i,1) + (stepa * sin(2*pi*(j/h)));
                            X2new = X(i,2) + (stepa * sin(2*pi*(j/h)));
                            if easom(X1new,X2new)<fitnessnew
                                x1new = X1new; x2new = X2new; fitnessnew = easom(X1new,X2new);
                                movefl = 1;
                            end
                        end
                        if movefl
                            X(i:1) = x1new; X(i:2) = x2new;
                            fitness(i) = fitnessnew;
                            m = X(:,1); n = X(:,2);
                            refreshdata
                            drawnow
                        end
                        if easom(X(i,1),X(i,2))<lead
                            lead = easom(X(i,1),X(i,2));
                            idx_lead = i;
                            brkfl=0;
                            break
                        end
                    end
                end
                T = T + 1;
            end
            step=3;
        elseif step==3
            % Calling
            distfl = 1;
            while distfl
                moveable_wolves = [];
                for i=1:N
                    dist(i) = sqrt((X(idx_lead,1)-X(i,1))^2 + (X(idx_lead,2)-X(i,2))^2);
                    if dist(i)>Lnear
                        moveable_wolves = [moveable_wolves i];
                    end
                end
                if isempty(moveable_wolves)
                    disfl=0;
                    break
                end
                for i=1:length(moveable_wolves)
                    x1old = X(moveable_wolves(i),1); x2old = X(moveable_wolves(i),2);
                    x1new = X(idx_lead,1); x2new = X(idx_lead,2);
                    X(moveable_wolves(i),1) = x1old + stepb * ((x1new-x1old)/abs(x1new-x1old));
                    X(moveable_wolves(i),2) = x2old + stepb * ((x2new-x2old)/abs(x2new-x2old));
                    m = X(:,1); n = X(:,2);
                    refreshdata
                    drawnow
                    if easom(X(moveable_wolves(i),1), X(moveable_wolves(i),2))<lead
                        lead = easom(X(moveable_wolves(i),1), X(moveable_wolves(i),2));
                        idx_lead = moveable_wolves(i);
                        step=2;
                        distfl = 0;
                        break
                    end
                end
            end
            if step==3
                step=4;
            end
        elseif step==4
            % Besieging
            for i=1:N
                if i~=idx_lead
                    x1old = X(i,1); x2old = X(i,2);
                    x1lead = X(idx_lead,1); x2lead = X(idx_lead,2);
                    x1new = x1old + (2*rand-1) * stepc * (abs(x1lead-x1old));
                    x2new = x2old + (2*rand-1) * stepc * (abs(x2lead-x2old));
                    if easom(x1new, x2new)<easom(x1old, x2old)
                        X(i,1) = x1new; X(i,2) = x2new;
                        m = X(:,1); n = X(:,2);
                        refreshdata
                        drawnow
                        if easom(x1new, x2new)<lead
                            lead = easom(x1new, x2new);
                            idx_lead = i;
                        end
                    end
                end
            end
            step=5;
        elseif step==5
            % stronger surviving renewing
            fitness = easom(X(:,1),X(:,2));
            [~, idx_sorted] = sort(fitness, "descend");
            R = randi([Rmin, Rmax], 1);
            for i=1:R
                    X(idx_sorted(i),1) = X(idx_lead, 1) * (1+(0.2*rand-0.1));
                    X(idx_sorted(i),2) = X(idx_lead, 2) * (1+(0.2*rand-0.1));
                    m = X(:,1); n = X(:,2);
                    refreshdata
                    drawnow
            end
            overfl=1;
        end
    end

    k = k + 1;
    fprintf("iteration: %d\n", k)
end

%%% Plotting fitness over time
figure;
plot(best);
hold on
plot(worst); plot(avg);
hold off
ylim([-1, 0])
title("Best, worst and average fitness");
xlabel("iterations","FontWeight", "bold"); ylabel("Fitness", "FontWeight", "bold");
legend("Best", "Worst", "Average", "Location", "best");