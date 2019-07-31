hold off
clear
clc
clf
%clears everything
hold on
width = 20;
height = 20; %define area
tEnd = 100; %time to run
I = 1; %# of infected
N = 200; %# number of agents
R = 0; % number of recovered
xs = zeros(N,1); %x positions
ys = zeros(N,1); %y positions
ixs = ones(N,1)*-1; %infected x positions
iys = ones(N,1)*-1; %infected y positions
rxs = ones(N,1)*-1; 
rys = ones(N,1)*-1; 
pM = 1; %probability of movement
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

tic


for count = 1:N %initializes the structure definitions of each agent
    population(count).name = count; %gives a label
    population(count).x = round(rand(1,1)*width); %defines x positions randomly
    population(count).y = round(rand(1,1)*height); %defines y position randomly
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
    clf %clears the figure
    hold on
    plot(xs,ys,'.','markersize',20) %plots health
    ixs(1:I) = xs(1:I); %defines sick xs
    iys(1:I) = ys(1:I); %defines sick ys
    rxs(1:R) = xs((N-R+1):N); %defines sick xs
    rys(1:R) = ys((N-R+1):N); %defines sick ys
    xlim([0,width])
    ylim([0,height])
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
    xlim([0,width])
    ylim([0,height])
    pause(0.001); %waits for display sake
    sick = [sick ; sum(ixs>0)]; %counts number of sick
    dead = [dead ; sum(ixs == -2)+sum(ixs == -5)];
    recovered = [recovered ; sum(ixs == -3)];
    hospitalized = [hospitalized ; sum(ixs == -4)];
    funeral = [funeral ; sum(ixs == -5)];
    healthy = [healthy ; N-sick(t)-dead(t)+recovered(t)-hospitalized(t)]; %counts number of healthy
end

toc
clf
plot(1:tEnd,healthy,1:tEnd,sick,1:tEnd,dead,1:tEnd,recovered,1:tEnd,hospitalized,1:tEnd,funeral)
legend('Susceptible','Infected','Dead','Recovered','Hospitalized','Funeral')

%% Heatmap 
hold off
clf
%clears everything
hold on

width = 20:20:200; %define area
pM = 0.1:0.1:1; %probability of movement
tEnd = 10;

I = 1; %# of infected
N = 200; %# number of agents
R = 0; % number of recovered
xs = zeros(N,1); %x positions
ys = zeros(N,1); %y positions
ixs = ones(N,1)*-1; %infected x positions
iys = ones(N,1)*-1; %infected y positions
rxs = ones(N,1)*-1; 
rys = ones(N,1)*-1;  %establish tracking

tic

for i = 1:length(width)
    for j = 1:length(pM)
        sick = [];
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
        ymtx(i,j,:,:) = [sick];
    end
end

toc