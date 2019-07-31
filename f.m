function I_real = f(constant,ttt)
   I_real = (constant(1)./((1+constant(2)).^ttt)).^ttt;
end 