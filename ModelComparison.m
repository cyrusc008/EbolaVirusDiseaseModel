%% Model Comparison for Compartmental Models
%Liberia = [3.4167e+03, 1.1236e+03, 3.4167e+03,1.1236e+03, 315.6397];
%Guinea  = [1.8204e+03,  369.5490, 1.8204e+03,369.5490, 170.7204];
%SierraL = [3.7857e+03, 1.7839e+03, 3.7857e+03,1.7839e+03, 539.3450];

Liberia = [ 3550, 1615, 3567, 1619, 418];
Guinea  = [ 6129, 7179, 5142, 4179, 2534];
SierraL = [ 1883, 547, 1900, 547, 272];

N_Lib = 7347;
N_Gui = 9796;
N_Sil = 2756;


PerLib = 100-(Liberia./N_Lib*100);
PerGui = 100-(Guinea./N_Gui*100);
PerSil = 100-(SierraL./N_Sil*100);

c = categorical({'a','b', 'c'});
y = [PerLib(5) PerGui(5) PerSil(5); PerLib(4) PerGui(4) PerSil(4); PerLib(1) PerGui(1) PerSil(1)];


figure;
hold on;
box on;
bar(c,y);
legend('Liberia 2014-16', 'Guinea 2014-16', 'Sierra Leone 2014-16','Location','northeast');
xlabel('Compartmental Models','FontSize',20);
ylabel('% Accuracy','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);
ylim([0 120]);

