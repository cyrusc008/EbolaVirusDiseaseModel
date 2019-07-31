%% BMED 4813 BHI: SIR MODEL V2(CYRUS) 
% 2014 Ebola outbreak in Liberia 
clear all, clc

%Load real datasets for total cases and deaths from CDC 
[filename, filepath, ~] = uigetfile('*.xlsx');
[~, headers, ~] = xlsread([filepath, filename], 1, '1:1');

tt = transpose(xlsread([filepath, filename], 1, 'A:A')); %day 0 is 25/3/2014
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

%Parameterization of the model 
N    = 3700;          %fit is 2000, normal is 1500000
E    = 32;
S    = N-E;

beta_IR = 0.250;         %fit is 0.170, normal is 0.160
beta_ID = 0.250;         %fit is 0.170, normal is 0.160
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


%Solving the ODEs
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

%Plot the SIR curves 
figure;
hold on;
box on;
plot(t-100,xa(:,1));
plot(t-100,xa(:,2));
plot(t-100,xa(:,3)+xa(:,4));
plot(t-100,xa(:,5)+xa(:,6));
plot(t-100,xa(:,7));
plot(t-100,xa(:,8));
plot(t-100,xa(:,9));
legend('Susceptibles','Exposed','Infectious','Hospitalized','Recovered','Funeral','Deceased');
xlabel('Time (Days)','FontSize',20);
ylabel('Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

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

%RMSE calculation 
ispan = 1:1500;
N = zeros(length(ispan),1);

for N = 1:15
    [t,xa]=ode45(f,tspan, [S E 0 0 0 0 0 0 0]);
    figure;
    plot(tt,xa(:,1));
%    I_eqn(:,i+1)  = 0*I_eqn(:,i)+cumsum(xa(:,3)+xa(:,4));
end 
        
    F=ode45(f,tspan, [S E 0 0 0 0 0 0 0]);
    I_eqn(i+1)  = 0*I_eqn(i) + cumsum(xa(:,3)+xa(:,4));
    RMSE(i+1) = 0*err(i)+(immse(I_real, transpose(I_eqn)))^0.5;