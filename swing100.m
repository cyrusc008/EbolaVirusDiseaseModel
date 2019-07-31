   function dydt= swing100(t,y)
   global Tm Te D H wB ;
   dydt=[wB.*y(2);(0.5/H).*(Tm-Te-(D.*y(2)))];