function y = QDEfitnessvalue1entropy_sir(x, table)
% This function is used to find the fitness value of each chromosome

finalim = x;
k = table(1,1);
l = table(1,2);
rmid = 90;
r1 = 5;  % Overwrites table(1,3), ensuring r1 is always 5
r2 = 5;  % Overwrites table(1,4), ensuring r2 is always 5

h = (rmid + r1) * 2 + 1;
w = (rmid + r2) * 2 + 1;

R = finalim(l - rmid - r1 : l + rmid + r1, k - rmid - r2 : k + rmid + r2);

y = entropy(R);

end