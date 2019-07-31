%% BMED 4813 BHI: SIR INITIAL MODEL (CYRUS) 
% 2014 Ebola outbreak in Liberia 

%Parameterization of the model 
N    = 4294000;          %CDC
I    = 846;              %WHO
R    = 735;              %WHO
S    = N-I+R;            %mortzlity rate
T    = 12;               %Rivers et al., 2014
M    = 0.5;              %Rivers et al., 2014

r_R  = 1/T;              %Recovery coefficient 
r_I  = M/S;              %Infectious coefficient 

%Solving the ODEs
f = @(t,x) [-r_I*x(1)*x(2);r_I*x(1)*x(2)-r_R*x(2);r_R*x(2)];
[t,xa]=ode45(f,[0 60], [S I R]);

%Plot the SIR curves 
figure;
hold on;
box on;
plot(t,xa(:,1));
plot(t,xa(:,2),'k');
plot(t,xa(:,3),'r');
legend('Susceptible','Infected','Recovered');
xlabel('Time (Days)','FontSize',20);
ylabel('Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

%Compare to real data (from WHO reports) 
tt      = [7 10 14 18 20 24 26 28 33 35 40 42 47 49 54 57];                
I_real  = [1871 2046 2081 2407 2710 3022 3280 3458 3696 3834 3924 4076 4249 4262 4665 4665];
I_eqn   = cumsum(xa(:,2));

figure;
hold on;
box on;
plot(t,I_eqn);
scatter(tt,I_real);
legend('Infected (SIR model)','Infected (real)');
xlabel('Time (Days)','FontSize',20);
ylabel('Cumulated Population','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);

%Nonlinear regression 
tt      = [7 10 14 18 20 24 26 28 33 35 40 42 47 49 54 57];                
I_real  = [1871 2046 2081 2407 2710 3022 3280 3458 3696 3834 3924 4076 4249 4262 4665 4665];
I_eqn   = cumsum(xa(:,2));

g = @(b,x) [x(1)*x(2)*(M/12)-b*x(2)-x(2)*M]; 
beta0 = [S I R];
mdl = fitnlm(transpose(tt),transpose(I_real),g,beta0); 