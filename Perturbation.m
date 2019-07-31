%% Model Perturbation 

factor         = [1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3];
beta_IR_ID_Lib = [1 1.0004 1.0463 1.1108 1.121 0.981573099 0.645933273 0.289109524];
beta_IR_ID_Gui = [1 0.996274885 0.997290178 0.979611523 0.909341944 0.758570315 0.541769949 0.33597077];
beta_IR_ID_SiL = [1 1.050536609 1.142913884 1.227578377 1.221418046 1.035268423 0.641345776 0.264345039];

beta_IR_avg    = [1 1.015737165 1.062168021 1.105996633 1.083919997 0.925137279 0.609682999 0.296475111];

%beta_F_Lib     = [1.005350885 1.000001822 0.996252716 0.994084319 0.993446295 0.994278245 0.9964937 1];
%beta_F_Gui     = [0.998134515 0.997072968 0.996280847 0.995888235 0.995996212 0.996683714 0.998006494 1];
%beta_F_SiL     = [1.065768148 1.051905952 1.039422627 1.02842082 1.018962889 1.011076996 1.004762315 1];

factor_2       = [1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7];

pie_Lib        = [1 0.955273958 0.915612178 0.880352268 0.848901602 0.820728686 0.795364147 0.772391528];
pie_Gui        = [1 0.956566448 0.916931011 0.880574979 0.847063995 0.816025695 0.787145646 0.760156514];
pie_SiL        = [1 0.964324497 0.932990466 0.905301271 0.880641045 0.858474586 0.838332412 0.819824767];

pie_avg        = [1 0.958721634 0.921844552 0.88874284 0.858868881 0.831742989 0.806947402 0.784124269];

e2_k2_Lib      = [1 0.967214523 0.938222919 0.91245443 0.88944184 0.868796739 0.85019412 0.833360965];
e2_k2_Gui      = [1 0.968076088 0.939518432 0.913828317 0.890591971 0.869477066 0.850205894 0.832546795];
e2_k2_SiL      = [1 0.972382525 0.948133348 0.926714772 0.907679633 0.890672427 0.875394805 0.86160635];

e2_k2_avg      = [1 0.969224379 0.941958233 0.91766584 0.895904481 0.876315411 0.858598273 0.842504703];

%gamma_Lib      = [1 0.99661856 0.994518691 0.993304875 0.992705454 0.992533059 0.992660912 0.992997932];
%gamma_Gui      = [1 0.997793994 0.99635852 0.995435449 0.994856669 0.99451251 0.994330038 0.994258801];
%gamma_SiL      = [1 1.003662908 1.007659754 1.011742372 1.015776913 1.019682966 1.02341878 1.026965309];

figure;
hold on;
box on;
%plot(factor,beta_IR_ID_Lib,factor,beta_IR_ID_Gui,factor,beta_IR_ID_SiL,factor_2,pie_Lib,factor_2,pie_Gui,factor_2,pie_SiL,factor_2,e2_k2_Lib,factor_2,e2_k2_Gui,factor_2,e2_k2_SiL,'LineWidth',3);
%plot(factor,beta_IR_ID_Lib,factor,beta_IR_ID_Gui,factor,beta_IR_ID_SiL,'LineWidth',3);
%plot(factor_2,pie_Lib,factor_2,pie_Gui,factor_2,pie_SiL,'LineWidth',3);
%plot(factor_2,e2_k2_Lib,factor_2,e2_k2_Gui,factor_2,e2_k2_SiL,'LineWidth',3);
%legend('Community Contact Rate, Liberia','Community Contact Rate, Guinea','Community Contact Rate, Sierra Leone','P of Seeking Hospitalization, Liberia','P of Seeking Hospitalization, Guinea','P of Seeking Hospitalization, Sierra Leone','Time Until Hospitalization, Liberia','Time Until Hospitalization, Guinea','Time Until Hospitalization, Sierra Leone');
scatter(factor,beta_IR_avg,'filled');
scatter(factor_2,pie_avg,'filled');
scatter(factor_2,e2_k2_avg,'filled');
legend('Community Contact Rate','P. of Seeking Hospitalization','Duration Until Hospitalization','location','southeast');
xlabel('Factor change in parameter','FontSize',20);
xlim([0 2.2]);
ylabel('Normalized total number of cases','FontSize',20);
set(gca, 'LineWidth',2,'FontSize',15);