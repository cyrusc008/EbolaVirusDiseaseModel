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

%Use the IDEA function to parameterize the SIR model 
ttt = linspace(0,6,length(tt));
constant = lsqcurvefit(@f, [0;0], ttt, I_real);

R_0 = constant(1);
d   = constant(2);

xfit = linspace(1,254,length(tt));
yfit = f(constant, xfit);   
