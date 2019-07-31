%% BMED 4813 BHI: SEIHRFD MODEL V2(CYRUS) 
% 2014 Ebola outbreak in Liberia 
clear all, clc

%% Load real datasets for total cases and deaths from CDC 
[filename, filepath, ~] = uigetfile('*.xlsx');
[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

tt = transpose(xlsread([filepath, filename], 1, 'A:A'));%day 0 is 25/3/2014
I_real = transpose(xlsread([filepath, filename], 1, 'E:E'));
D_real = transpose(xlsread([filepath, filename], 1, 'F:F'));

figure;
hold on;
box on;
scatter(tt,I_real);
scatter(tt,D_real);
legend('Infected Compartment (real)', 'Dead Compartment (real)');
xlabel('Time (Days)','FontSize',20);
ylabel('Number of People','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

%% Parameterization of the model 
N    = 7995;          %fit is 3900, normal is 1500000
E    = 32;
I    = 6;
H    = 1.4;
R    = 2.4;
F    = 3.6;
D    = 5;
S    = N-E-2*I-2*H-R-F-D;


beta_IR = 0.160510836;         %fit is 0.230, normal is 0.160
beta_ID = 0.161095084;         %fit is 0.230, normal is 0.160
beta_HR = 0.067663422;
beta_HD = 0.064597294;
beta_F  = 0.499593564;
theta   = 0.480964596;
alpha   = 0.079288266;
e_1     = 0.059975698;
e_2     = 0.274125202;
k_2     = 0.287440036;
k_1     = 0.071184955;
pie     = 0.138445185;
roe     = 0.067340937;
delta   = 0.105670587;
gamma   = 0.534135304;

%% Solving the ODEs
tspan = tt;
f = @(t,x) [-(1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8));
            (1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8))-alpha*x(2);
            (1-theta)*alpha*x(2)-(1-pie)*e_1*x(3)-pie*e_2*x(3);
            theta*alpha*x(2)-(1-pie)*k_1*x(4)-pie*k_2*x(4);
            pie*e_2*x(3)-roe*x(5);
            pie*k_2*x(4)-delta*x(6);
            (1-pie)*e_1*x(3)+roe*x(5);
            (1-pie)*k_1*x(4)-gamma*x(8);
            gamma*x(8)+delta*x(6)];
[t,xa]=ode45(f,tspan, [S E I I H H R F D]);

%% Plot the SIR curves 
figure;
hold on;
box on;
plot(t,xa(:,1),'LineWidth',3);
plot(t,xa(:,2),'LineWidth',3);
plot(t,xa(:,3)+xa(:,4),'LineWidth',3);
plot(t,xa(:,5)+xa(:,6),'LineWidth',3);
plot(t,xa(:,7),'LineWidth',3);
plot(t,xa(:,8),'LineWidth',3);
plot(t,xa(:,9),'LineWidth',3);
legend('Susceptibles','Exposed','Infectious','Hospitalized','Recovered','Funeral','Deceased');
xlabel('Time (Days)','FontSize',20);
ylabel('Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

I_eqn   = cumsum(xa(:,3)+xa(:,4));
RMSE = (immse(I_real, transpose(I_eqn)))^0.5

figure;
hold on;
box on;
plot(linspace(min(tt),max(tt),length(I_eqn)),I_eqn,'LineWidth',3);
scatter(tt,I_real);
legend('Infected (SIR model)','Infected (real)');
xlabel('Time (Days)','FontSize',20);
ylabel('Cumulative Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);
xlim([0 700]);
ylim([0 20000]);

%% RMSE minimization for N and contact rate
N = linspace(4000,10000,16);
beta_IR = linspace(0.1,0.5,16);
tic

for i = 1:length(N)
            for j = 1:length(beta_IR)
                S(i) = N(i)-E-2*I-2*H-R-F-D;
                f = @(t,x) [-(1./N(i)).*(beta_IR(j).*x(1).*x(3)+beta_IR(j).*x(1).*x(4)+beta_HR.*x(1).*x(5)+beta_HD.*x(1).*x(6)+beta_F.*x(1).*x(8));
                            (1./N(i)).*(beta_IR(j).*x(1).*x(3)+beta_IR(j).*x(1).*x(4)+beta_HR.*x(1).*x(5)+beta_HD.*x(1).*x(6)+beta_F.*x(1).*x(8))-alpha.*x(2);
                            (1-theta).*alpha.*x(2)-(1-pie).*e_1.*x(3)-pie.*e_2.*x(3);
                            theta.*alpha.*x(2)-(1-pie).*k_1.*x(4)-pie.*k_2.*x(4);
                            pie.*e_2.*x(3)-roe.*x(5);
                            pie.*k_2.*x(4)-delta.*x(6);
                            (1-pie).*e_1.*x(3)+roe.*x(5);
                            (1-pie).*k_1.*x(4)-gamma.*x(8);
                            gamma.*x(8)+delta.*x(6)];
                [t,xa]=ode45(f,tspan, [S(i) E I I H H R F D]);
                ymtx(i,j,:,:) = xa;
            end 
end

toc

for k = 1:length(N)
    for l = 1:length(beta_IR)
        xaa   = squeeze(ymtx(k,l,:,:));
        I_eqn = cumsum(xaa(:,3)+xaa(:,4)); 
        RMSE(k,l) = abs((N(k)-((immse(I_real, transpose(I_eqn)))^0.5))/N(k)*100);
    end
end 

figure; 
h = heatmap(round(transpose(beta_IR),2),round(transpose(N),0),RMSE);
h.Title = 'Sierra Leone 2014-16';
h.XLabel = 'Community Contact Rate';
h.YLabel = 'Total Population';
 h.Colormap = autumn(100);
h.FontSize = 20
h.GridVisible = 'off'
h.ColorLimits = [0 100];
%h.ColorLimits = [find(RMSE == min(RMSE(:))) find(RMSE == max(RMSE(:)))];

%% RMSE minimization for seeking hos and time until hospitalization 
pie = linspace(0.1,0.5,16);
e_2 = linspace(0.1,0.5,16);
tic

for i = 1:length(pie)
            for j = 1:length(e_2)
                S = N-E-2*I-2*H-R-F-D;
                f = @(t,x) [-(1./N).*(beta_IR.*x(1).*x(3)+beta_IR.*x(1).*x(4)+beta_HR.*x(1).*x(5)+beta_HD.*x(1).*x(6)+beta_F.*x(1).*x(8));
                            (1./N).*(beta_IR.*x(1).*x(3)+beta_IR.*x(1).*x(4)+beta_HR.*x(1).*x(5)+beta_HD.*x(1).*x(6)+beta_F.*x(1).*x(8))-alpha.*x(2);
                            (1-theta).*alpha.*x(2)-(1-pie(i)).*e_1.*x(3)-pie(i).*e_2(j).*x(3);
                            theta.*alpha.*x(2)-(1-pie(i)).*k_1.*x(4)-pie(i).*e_2(j).*x(4);
                            pie(i).*e_2(j).*x(3)-roe.*x(5);
                            pie(i).*e_2(j).*x(4)-delta.*x(6);
                            (1-pie(i)).*e_1.*x(3)+roe.*x(5);
                            (1-pie(i)).*k_1.*x(4)-gamma.*x(8);
                            gamma.*x(8)+delta.*x(6)];
                [t,xa]=ode45(f,tspan, [S E I I H H R F D]);
                ymtx(i,j,:,:) = xa;
            end 
end

toc

for k = 1:length(pie)
    for l = 1:length(e_2)
        xaa   = squeeze(ymtx(k,l,:,:));
        I_eqn = cumsum(xaa(:,3)+xaa(:,4)); 
        RMSE(k,l) = abs((N-((immse(I_real, transpose(I_eqn)))^0.5))/N*100);
    end
end 

figure; 
h = heatmap(round(transpose(e_2),2),round(transpose(pie),2),RMSE);
h.Title = 'Sierra Leone 2014-16';
h.XLabel = 'Time Until Hospitalization';
h.YLabel = 'P. of Seeking Hospitals';
 h.Colormap = autumn(100);
h.FontSize = 20
h.GridVisible = 'off'
h.ColorLimits = [0 100];
%h.ColorLimits = [find(RMSE == min(RMSE(:))) find(RMSE == max(RMSE(:)))];

%%
beta_F = linspace(0.14925372,0.489,8);

tic

for i = 1:length(beta_F)
                S = N-E-2*I-2*H-R-F-D;
                f = @(t,x) [-(1./N).*(beta_IR.*x(1).*x(3)+beta_IR.*x(1).*x(4)+beta_HR.*x(1).*x(5)+beta_HD.*x(1).*x(6)+beta_F(i).*x(1).*x(8));
                            (1./N).*(beta_IR.*x(1).*x(3)+beta_IR.*x(1).*x(4)+beta_HR.*x(1).*x(5)+beta_HD.*x(1).*x(6)+beta_F(i).*x(1).*x(8))-alpha.*x(2);
                            (1-theta).*alpha.*x(2)-(1-pie).*e_1.*x(3)-pie.*e_2.*x(3);
                            theta.*alpha.*x(2)-(1-pie).*k_1.*x(4)-pie.*k_2.*x(4);
                            pie.*e_2.*x(3)-roe.*x(5);
                            pie.*k_2.*x(4)-delta.*x(6);
                            (1-pie).*e_1.*x(3)+roe.*x(5);
                            (1-pie).*k_1.*x(4)-gamma.*x(8);
                            gamma.*x(8)+delta.*x(6)];
                [t,xa]=ode45(f,tspan, [S E I I H H R F D]);
                ymtx(i,:,:) = xa;
end

toc

for k = 1:length(beta_F)
        xaa     = squeeze(ymtx(k,:,:));
        I_eqn   = cumsum(xaa(:,3)+xaa(:,4)); 
        Max(k)  = max(I_eqn);
        RMSE(k) = (immse(I_real, transpose(I_eqn)))^0.5;
end 

Maxx = Max/Max(length(beta_F)); 