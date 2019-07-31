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
N    = 3900;          %fit is 3900, normal is 1500000
E    = 32;
S    = N-E;

beta_IR = 0.230;         %fit is 0.230, normal is 0.160
beta_ID = 0.230;         %fit is 0.230, normal is 0.160
beta_HR = 0.062;
beta_HD = 0.062;
beta_F  = 0.489;
theta   = 0.45049;
alpha   = 0.088883;
e_1     = 0.066667;
e_2     = 0.3086;
k_2     = 0.3086;
k_1     = 0.07513148;
pie     = 0.197;
roe     = 0.06297229;
delta   = 0.09930487;
gamma   = 0.4975124;


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
[t,xa]=ode45(f,tspan, [S E 0 0 0 0 0 0 0]);

%% Plot the SIR curves 
figure;
hold on;
box on;
plot(t-100,xa(:,1),'LineWidth',3);
plot(t-100,xa(:,2),'LineWidth',3);
plot(t-100,xa(:,3)+xa(:,4),'LineWidth',3);
plot(t-100,xa(:,5)+xa(:,6),'LineWidth',3);
plot(t-100,xa(:,7),'LineWidth',3);
plot(t-100,xa(:,8),'LineWidth',3);
plot(t-100,xa(:,9),'LineWidth',3);
legend('Susceptibles','Exposed','Infectious','Hospitalized','Recovered','Funeral','Deceased');
xlabel('Time (Days)','FontSize',20);
ylabel('Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

I_eqn   = cumsum(xa(:,3)+xa(:,4));

figure;
hold on;
box on;
plot(t,I_eqn,'LineWidth',3);
scatter(tt,I_real);
legend('Infected (SIR model)','Infected (real)');
xlabel('Time (Days)','FontSize',20);
ylabel('Cumulative Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);
xlim([0 260]);

%% RMSE minimization for N 
N = 1500;
S = N-E;
f = @(t,x) [-(1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8));
            (1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8))-alpha*x(2);
            (1-theta)*alpha*x(2)-(1-pie)*e_1*x(3)-pie*e_2*x(3);
            theta*alpha*x(2)-(1-pie)*k_1*x(4)-pie*k_2*x(4);
            pie*e_2*x(3)-roe*x(5);
            pie*k_2*x(4)-delta*x(6);
            (1-pie)*e_1*x(3)+roe*x(5);
            (1-pie)*k_1*x(4)-gamma*x(8);
            gamma*x(8)+delta*x(6)];
[t,xa]=ode45(f,tspan, [S E 0 0 0 0 0 0 0]);
I_eqn   = cumsum(xa(:,3)+xa(:,4));
RMSE = (immse(I_real, transpose(I_eqn)))^0.5;
disp(RMSE);

%[filename, filepath, ~] = uigetfile('*.xlsx');
%[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

%NN = transpose(xlsread([filepath, filename], 1, 'A:A')); 
%N_rmse = transpose(xlsread([filepath, filename], 1, 'B:B')); 

%figure;
%hold on;
%box on;
%plot(NN,N_rmse,'LineWidth',3);
%xlabel('N','FontSize',20);
%ylabel('RMSE','FontSize',20);
%set(gca, 'xscale','log','yscale','log','LineWidth',2,'FontSize',15);

%% RMSE minimization for beta IR/ID
beta_IR = 0.17;         
beta_ID = 0.17;         

f = @(t,x) [-(1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8));
            (1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8))-alpha*x(2);
            (1-theta)*alpha*x(2)-(1-pie)*e_1*x(3)-pie*e_2*x(3);
            theta*alpha*x(2)-(1-pie)*k_1*x(4)-pie*k_2*x(4);
            pie*e_2*x(3)-roe*x(5);
            pie*k_2*x(4)-delta*x(6);
            (1-pie)*e_1*x(3)+roe*x(5);
            (1-pie)*k_1*x(4)-gamma*x(8);
            gamma*x(8)+delta*x(6)];
[t,xa]=ode45(f,tspan, [S E 0 0 0 0 0 0 0]);
I_eqn   = cumsum(xa(:,3)+xa(:,4));
RMSE = (immse(I_real, transpose(I_eqn)))^0.5;
disp(RMSE);

[filename, filepath, ~] = uigetfile('*.xlsx');
[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

beta_IRID = transpose(xlsread([filepath, filename], 1, 'E:E')); 
beta_rmse = transpose(xlsread([filepath, filename], 1, 'F:F')); 

figure;
hold on;
box on;
plot(beta_IRID,beta_rmse,'LineWidth',3);
xlabel('beta_IR/ID','FontSize',20);
ylabel('RMSE','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

%% RMSE minimization for overall SEIHRFD Model 
N = 14000;
beta_IR = 0.40;
beta_ID = 0.40;   
S    = N-E;

f = @(t,x) [-(1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8));
            (1/N)*(beta_IR*x(1)*x(3)+beta_ID*x(1)*x(4)+beta_HR*x(1)*x(5)+beta_HD*x(1)*x(6)+beta_F*x(1)*x(8))-alpha*x(2);
            (1-theta)*alpha*x(2)-(1-pie)*e_1*x(3)-pie*e_2*x(3);
            theta*alpha*x(2)-(1-pie)*k_1*x(4)-pie*k_2*x(4);
            pie*e_2*x(3)-roe*x(5);
            pie*k_2*x(4)-delta*x(6);
            (1-pie)*e_1*x(3)+roe*x(5);
            (1-pie)*k_1*x(4)-gamma*x(8);
            gamma*x(8)+delta*x(6)];
[t,xa]=ode45(f,tspan, [S E 0 0 0 0 0 0 0]);
I_eqn   = cumsum(xa(:,3)+xa(:,4));

figure;
hold on;
box on;
plot(t,I_eqn);
scatter(tt,I_real);
legend('Infected (SIR model)','Infected (real)');
xlabel('Time (Days)','FontSize',20);
ylabel('Cumulative Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);
xlim([0 260]);
RMSE = (immse(I_real, transpose(I_eqn)))^0.5;
disp(RMSE);

[filename, filepath, ~] = uigetfile('*.xlsx');
[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

NNN = transpose(xlsread([filepath, filename], 1, 'G:G')); 
bb_IRID = transpose(xlsread([filepath, filename], 1, 'H:H')); 
rmse_T = transpose(xlsread([filepath, filename], 1, 'I:I')); 

figure;
hold on;
box on;
scatter3(NNN,bb_IRID,rmse_T,'filled');
xlabel('N','FontSize',20);
ylabel('beta IR/ID','FontSize',20);
zlabel('RMSE', 'FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);
view(40,35)
