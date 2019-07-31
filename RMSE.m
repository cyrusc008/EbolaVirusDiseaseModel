function err = RMSE(tspan,I_real_2,par, S, E, I, H, R, F, D)
[t,xa] = ode45(@(t,X)SIRSolver(t,X,par), tspan, [S E I I H H R F D]); 
err = sum(sum((I_real_2-transpose(cumsum(xa(:,3)+xa(:,4)))).^2));
end