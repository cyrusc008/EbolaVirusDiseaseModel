hold off
clear
clc
clf
%clears everything
hold on
width = 155;
height = 155; %define area
tEnd = 100; %time to run
I = 1; %# of infected
N = 100; %# number of agents
R = 0; % number of recovered
xs = zeros(N,1); %x positions
ys = zeros(N,1); %y positions
ixs = ones(N,1)*-1; %infected x positions
iys = ones(N,1)*-1; %infected y positions
rxs = ones(N,1)*-1; 
rys = ones(N,1)*-1; 
pM = 0.55; %probability of movement
pT = 0.7; %probability of transmission
pDU = 0.70; %probability of death (unhospitalized)
pDH = 0.45; %probability of death (hospitalized) 
pH = 0.2 ; %probability of seeking hospitalization 
pF = (1-pH)*pDU ; %probability of funeral, and by extension death
sick = [];
dead = [];
recovered = [];
hospitalized = [];
funeral = [];
healthy = []; %establish tracking
inter = [];

[filename, filepath, ~] = uigetfile('*.xlsx');
[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

tt = transpose(xlsread([filepath, filename], 1, 'A:A'));%day 0 is 25/3/2014
I_real = transpose(xlsread([filepath, filename], 1, 'E:E'));

tic


for count = 1:N %initializes the structure definitions of each agent
    population(count).name = count; %gives a label
    population(count).x = round(rand(1,1)*width); %defines x positions randomly
    population(count).y = round(rand(1,1)*height); %defines y position randomly
    population(count).incubation = 0; %defines every agent as healthy
    xs(population(count).name) = population(count).x; %sets the xs value to the x position
    ys(population(count).name) = population(count).y; %sets the ys value to the y position
end
for count = 1:I
    population(count).incubation = floor(rand(1,1)*28+2); %sets first x population to be sick
end
ixs(1:I) = xs(1:I); %defines sick xs
iys(1:I) = ys(1:I); %defines sick ys
rxs(1:R) = xs((N-R+1):N); %defines sick xs
rys(1:R) = ys((N-R+1):N); %defines sick ys
for t = 1:tEnd %loops for time
    for count = 1:N %loops per agent
        m = rand(1,1); %chance to move
        if m<pM %determines if moves
            m = floor(rand(1,1)*8+1); %determines direction
            %m = 8;
            switch m %goes into case based on direction
                case 1
                    if population(count).x < width
                        population(count).x = population(count).x + 1;
                        xs(count) = population(count).x;
                    else
                        population(count).x = population(count).x - 1;
                        xs(count) = population(count).x;
                    end
                case 2 %example of down and to the right
                    if (population(count).x < width)&&(population(count).y > 0) %checks if not at either bound
                        population(count).x = population(count).x + 1; %moves the x
                        population(count).y = population(count).y - 1; %moves the y
                        xs(count) = population(count).x; %sets xs
                        ys(count) = population(count).y; %sets ys
                    elseif population(count).x < width %checks if it is at just the y bound
                        population(count).y = population(count).y + 1; %bounces the opposite direction
                        ys(count) = population(count).y; %sets ys
                    else
                        population(count).x = population(count).x - 1; %bounces x
                        xs(count) = population(count).x; %sets x
                    end
                case 3
                    if population(count).y > 0
                        population(count).y = population(count).y - 1;
                        ys(count) = population(count).y;
                    else
                        population(count).y = population(count).y + 1;
                        ys(count) = population(count).y;
                    end
                case 4
                    if (population(count).x > 0)&&(population(count).y > 0)
                        population(count).x = population(count).x - 1;
                        population(count).y = population(count).y - 1;
                        xs(count) = population(count).x;
                        ys(count) = population(count).y;
                    elseif population(count).x > 0
                        population(count).y = population(count).y + 1;
                        ys(count) = population(count).y;
                    else
                        population(count).x = population(count).x + 1;
                        xs(count) = population(count).x;
                    end
                case 5
                    if population(count).x > 0
                        population(count).x = population(count).x - 1;
                        xs(count) = population(count).x;
                    else
                        population(count).x = population(count).x + 1;
                        xs(count) = population(count).x;
                    end
                case 6
                    if (population(count).x > 0)&&(population(count).y < height)
                        population(count).x = population(count).x - 1;
                        population(count).y = population(count).y + 1;
                        xs(count) = population(count).x;
                        ys(count) = population(count).y;
                    elseif population(count).x > 0
                        population(count).y = population(count).y - 1;
                        ys(count) = population(count).y;
                    else
                        population(count).x = population(count).x + 1;
                        xs(count) = population(count).x;
                    end
                case 7
                    if population(count).y < height
                        population(count).y = population(count).y + 1;
                        ys(count) = population(count).y;
                    else
                        population(count).y = population(count).y - 1;
                        ys(count) = population(count).y;
                    end
                case 8
                    if (population(count).x < width)&&(population(count).y < height)
                        population(count).x = population(count).x + 1;
                        population(count).y = population(count).y + 1;
                        xs(count) = population(count).x;
                        ys(count) = population(count).y;
                    elseif population(count).x < width
                        population(count).y = population(count).y - 1;
                        ys(count) = population(count).y;
                    else
                        population(count).x = population(count).x - 1;
                        xs(count) = population(count).x;
                    end
            end
        end
    end
    ixs(1:I) = xs(1:I); %defines sick xs
    iys(1:I) = ys(1:I); %defines sick ys
    rxs(1:R) = xs((N-R+1):N); %defines sick xs
    rys(1:R) = ys((N-R+1):N); %defines sick ys
    for count = 1:N %loops per agent
        for pos = 1:N %loops per sick
            if (xs(count) == ixs(pos))&&(ys(count)==iys(pos))&&population(count).incubation <1 %checks if positions overlap
                T = rand(1,1); %random to pass on disease
                if T<pT %determines if disease passes
                    population(count).incubation = round(rand(1,1)*19+2); %sets incubation
                    rxs(1:sum(ixs == -3)) = xs((N-sum(ixs == -3)+1):N); %defines sick xs
                    rys(1:sum(ixs == -3)) = ys((N-sum(ixs == -3)+1):N); %defines sick ys
                end
            end
        end
           
    end
    for count = 1:N
       if population(count).incubation>0 %checks if infected
           ixs(count) = xs(count); %adds to infected pool
           iys(count) = ys(count);
       else
           ixs(count) = -1; %removes from infected pool if not infected
           iys(count) = -1;
       end
       for state = 1:N
           D = rand(1,1);
           if (D<pDU) && (population(count).incubation > 0) && (population(count).incubation < 5)
               ixs(count) = -2;
           elseif (D>pDU) && (population(count).incubation > 0) && (population(count).incubation < 5)
               rxs(count) = ixs(count); 
               rys(count) = iys(count); 
               ixs(count) = -3;
           end       
           if (D<pDH) && (population(count).incubation > 0) && (population(count).incubation < 5)
               ixs(count) = -2;
           elseif (D>pDH) && (population(count).incubation > 0) && (population(count).incubation < 5)
               rxs(count) = ixs(count); 
               rys(count) = iys(count); 
               ixs(count) = -3;
           end                
           H = rand(1,1);  
           if (H<pH) && (population(count).incubation > 15) && (population(count).incubation < 20)
               ixs(count) = -4;
           end
           F = rand(1,1);
           if (F<pF) && (population(count).incubation > 0) && (population(count).incubation < 5)
               ixs(count) = -5;
           end
       end
       population(count).incubation = population(count).incubation -1; %ticks down incubation
    end
    sick = [sick ; sum(ixs>0)]; %counts number of sick
    dead = [dead ; sum(ixs == -2)+sum(ixs == -5)];
    recovered = [recovered ; sum(ixs == -3)];
    hospitalized = [hospitalized ; sum(ixs == -4)];
    funeral = [funeral ; sum(ixs == -5)];
    healthy = [healthy ; N-sick(t)-dead(t)+recovered(t)-hospitalized(t)]; %counts number of healthy
end

toc
clf

for i = 1:length(tt)
    xaa(i) = transpose(sick(tt(i)));
end

I_eqn = cumsum(xaa);

plot(1:tEnd,healthy,1:tEnd,sick,1:tEnd,dead,1:tEnd,recovered,1:tEnd,hospitalized,1:tEnd,funeral)
legend('Susceptible','Infected','Dead','Recovered','Hospitalized','Funeral')

figure;
hold on;
box on;
plot(tt,I_eqn,'LineWidth',3);
scatter(tt,I_real);
legend('Infected (SIR model)','Infected (real)');
xlabel('Time (Days)','FontSize',20);
ylabel('Cumulative Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);
xlim([0 20]);

RMSE = (immse(I_real, I_eqn))^0.5;

%% Heatmap 
hold off
clf
%clears everything
hold on

[filename, filepath, ~] = uigetfile('*.xlsx');
[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

tt = transpose(xlsread([filepath, filename], 1, 'A:A'));%day 0 is 25/3/2014
I_real = transpose(xlsread([filepath, filename], 1, 'E:E'));

boundary_high_width = 2000;
boundary_low_width = 20;
stepsize = 100;
width = linspace(boundary_low_width,boundary_high_width,stepsize); %define area

boundary_high_pM = 1;
boundary_low_pM = 0.01;
pM = linspace(boundary_low_pM,boundary_high_pM,stepsize); %define area
tEnd = 100;

I = 1; %# of infected
N = 100; %# number of agents
R = 0; % number of recovered
xs = zeros(N,1); %x positions
ys = zeros(N,1); %y positions
ixs = ones(N,1)*-1; %infected x positions
iys = ones(N,1)*-1; %infected y positions
rxs = ones(N,1)*-1; 
rys = ones(N,1)*-1;  %establish tracking
ymtx = zeros(length(width),length(pM));
tic

for i = 1:length(width)
    inter = [];
    for j = 1:length(pM)
        sick = [];
        sick2 = [];
        dead = [];
        recovered = [];
        hospitalized = [];
        funeral = [];
        healthy = [];
        for count = 1:N %initializes the structure definitions of each agent
            population(count).name = count; %gives a label
            population(count).x = round(rand(1,1)*width(i)); %defines x positions randomly
            population(count).y = round(rand(1,1)*width(i)); %defines y position randomly
            population(count).incubation = 0; %defines every agent as healthy
            xs(population(count).name) = population(count).x; %sets the xs value to the x position
            ys(population(count).name) = population(count).y; %sets the ys value to the y position
        end
        plot(xs,ys,'.','markersize',20); %starts plot of everyone
        for count = 1:I
            population(count).incubation = floor(rand(1,1)*28+2); %sets first x population to be sick
        end
        ixs(1:I) = xs(1:I); %defines sick xs
        iys(1:I) = ys(1:I); %defines sick ys
        rxs(1:R) = xs((N-R+1):N); %defines sick xs
        rys(1:R) = ys((N-R+1):N); %defines sick ys
        plot(ixs(1:I)>0,iys(1:I)>0,'r.','markersize',20) %plots sick
        plot(rxs(1:R)>0,rys(1:R)>0,'g.','markersize',20) %plots recovered
        for t = 1:tEnd %loops for time
            for count = 1:N %loops per agent
                m = rand(1,1); %chance to move
                if m<pM(j) %determines if moves
                    m = floor(rand(1,1)*8+1); %determines direction
                    %m = 8;
                    switch m %goes into case based on direction
                        case 1
                            if population(count).x < width(i)
                                population(count).x = population(count).x + 1;
                                xs(count) = population(count).x;
                            else
                                population(count).x = population(count).x - 1;
                                xs(count) = population(count).x;
                            end
                        case 2 %example of down and to the right
                            if (population(count).x < width(i))&&(population(count).y > 0) %checks if not at either bound
                                population(count).x = population(count).x + 1; %moves the x
                                population(count).y = population(count).y - 1; %moves the y
                                xs(count) = population(count).x; %sets xs
                                ys(count) = population(count).y; %sets ys
                            elseif population(count).x < width(i) %checks if it is at just the y bound
                                population(count).y = population(count).y + 1; %bounces the opposite direction
                                ys(count) = population(count).y; %sets ys
                            else
                                population(count).x = population(count).x - 1; %bounces x
                                xs(count) = population(count).x; %sets x
                            end
                        case 3
                            if population(count).y > 0
                                population(count).y = population(count).y - 1;
                                ys(count) = population(count).y;
                            else
                                population(count).y = population(count).y + 1;
                                ys(count) = population(count).y;
                            end
                        case 4
                            if (population(count).x > 0)&&(population(count).y > 0)
                                population(count).x = population(count).x - 1;
                                population(count).y = population(count).y - 1;
                                xs(count) = population(count).x;
                                ys(count) = population(count).y;
                            elseif population(count).x > 0
                                population(count).y = population(count).y + 1;
                                ys(count) = population(count).y;
                            else
                                population(count).x = population(count).x + 1;
                                xs(count) = population(count).x;
                            end
                        case 5
                            if population(count).x > 0
                                population(count).x = population(count).x - 1;
                                xs(count) = population(count).x;
                            else
                                population(count).x = population(count).x + 1;
                                xs(count) = population(count).x;
                            end
                        case 6
                            if (population(count).x > 0)&&(population(count).y < width(i))
                                population(count).x = population(count).x - 1;
                                population(count).y = population(count).y + 1;
                                xs(count) = population(count).x;
                                ys(count) = population(count).y;
                            elseif population(count).x > 0
                                population(count).y = population(count).y - 1;
                                ys(count) = population(count).y;
                            else
                                population(count).x = population(count).x + 1;
                                xs(count) = population(count).x;
                            end
                        case 7
                            if population(count).y < width(i)
                                population(count).y = population(count).y + 1;
                                ys(count) = population(count).y;
                            else
                                population(count).y = population(count).y - 1;
                                ys(count) = population(count).y;
                            end
                        case 8
                            if (population(count).x < width(i))&&(population(count).y < width(i))
                                population(count).x = population(count).x + 1;
                                population(count).y = population(count).y + 1;
                                xs(count) = population(count).x;
                                ys(count) = population(count).y;
                            elseif population(count).x < width(i)
                                population(count).y = population(count).y - 1;
                                ys(count) = population(count).y;
                            else
                                population(count).x = population(count).x - 1;
                                xs(count) = population(count).x;
                            end
                    end
                end
            end
            clf %clears the figure
            hold on
            plot(xs,ys,'.','markersize',20) %plots health
            ixs(1:I) = xs(1:I); %defines sick xs
            iys(1:I) = ys(1:I); %defines sick ys
            rxs(1:R) = xs((N-R+1):N); %defines sick xs
            rys(1:R) = ys((N-R+1):N); %defines sick ys
            xlim([0,width(i)])
            ylim([0,width(i)])
            for count = 1:N %loops per agent
                for pos = 1:N %loops per sick
                    if (xs(count) == ixs(pos))&&(ys(count)==iys(pos))&&population(count).incubation <1 %checks if positions overlap
                        T = rand(1,1); %random to pass on disease
                        if T<pT %determines if disease passes
                            population(count).incubation = round(rand(1,1)*19+2); %sets incubation
                            rxs(1:sum(ixs == -3)) = xs((N-sum(ixs == -3)+1):N); %defines sick xs
                            rys(1:sum(ixs == -3)) = ys((N-sum(ixs == -3)+1):N); %defines sick ys
                        end
                    end
                end

            end
            for count = 1:N
               if population(count).incubation>0 %checks if infected
                   ixs(count) = xs(count); %adds to infected pool
                   iys(count) = ys(count);
               else
                   ixs(count) = -1; %removes from infected pool if not infected
                   iys(count) = -1;
               end
               for state = 1:N
                   D = rand(1,1);
                   if (D<pDU) && (population(count).incubation > 0) && (population(count).incubation < 5)
                       ixs(count) = -2;
                   elseif (D>pDU) && (population(count).incubation > 0) && (population(count).incubation < 5)
                       rxs(count) = ixs(count); 
                       rys(count) = iys(count); 
                       ixs(count) = -3;
                   end       
                   if (D<pDH) && (population(count).incubation > 0) && (population(count).incubation < 5)
                       ixs(count) = -2;
                   elseif (D>pDH) && (population(count).incubation > 0) && (population(count).incubation < 5)
                       rxs(count) = ixs(count); 
                       rys(count) = iys(count); 
                       ixs(count) = -3;
                   end                
                   H = rand(1,1);  
                   if (H<pH) && (population(count).incubation > 15) && (population(count).incubation < 20)
                       ixs(count) = -4;
                   end
                   F = rand(1,1);
                   if (F<pF) && (population(count).incubation > 0) && (population(count).incubation < 5)
                       ixs(count) = -5;
                   end
               end
               population(count).incubation = population(count).incubation -1; %ticks down incubation
            end
            plot(ixs(ixs>0),iys(ixs>0),'r.','markersize',20) %plots again
            plot(rxs(rxs>0),rys(rxs>0),'g.','markersize',20) %plots again
            xlim([0,width(i)])
            ylim([0,width(i)])
            sick = [sick ; sum(ixs>0)]; %counts number of sick
            dead = [dead ; sum(ixs == -2)+sum(ixs == -5)];
            recovered = [recovered ; sum(ixs == -3)];
            hospitalized = [hospitalized ; sum(ixs == -4)];
            funeral = [funeral ; sum(ixs == -5)];
            healthy = [healthy ; N-sick(t)-dead(t)+recovered(t)-hospitalized(t)]; %counts number of healthy
        end
        inter = [inter,sick];
    end
    ymtx = cat(3,ymtx,inter);
end

toc

for i = 1:length(tt)
    for j = 1:length(tt)
        for k = 1:length(tt)
            xaa(i,j,k) = ymtx(tt(i),tt(j),tt(k)) 
        end
    end
end

for m = 1:length(tt)
    for n = 1:length(tt)
        xaa_final   = squeeze(xaa(m,n,:,:));
        I_eqn = cumsum(xaa_final(:,1));
        RMSE(m,n) = (immse(I_real, transpose(I_eqn)))^0.5;
    end
end

figure; 
h = heatmap(linspace(boundary_low_width,boundary_high_width,length(tt)),linspace(boundary_low_pM,boundary_high_pM,length(tt)),RMSE);
h.Title = 'Heat-map for Parameter Selection: Liberia 2014';
h.XLabel = 'Width/Height';
h.YLabel = 'Probability of Movement';
h.Colormap = autumn(100);
h.FontSize = 15
h.GridVisible = 'off'
%h.ColorLimits = [find(RMSE == min(RMSE(:))) find(RMSE == max(RMSE(:)))];