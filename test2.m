x = [1 2 3 4 5 6 7 8 9 10];
y = [2 4 6 8 10 12 14 16 18 20];
a = [1 4 8];
b = [3 5 7];
z = zeros(1,length(b));

for i = 1:length(b)
    z(i) = y(a(i)); 
end