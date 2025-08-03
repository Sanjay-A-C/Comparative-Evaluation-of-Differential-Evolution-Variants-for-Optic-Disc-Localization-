function [x1 y1 width height] = squaretable(optimqtable1)
%This program is used to find the OD centre from the best chromosome
k=optimqtable1(1,1);
l=optimqtable1(1,2);
rmid=90;
r1=optimqtable1(1,3);
r2=optimqtable1(1,4);

h=(rmid+r1)*2+1;
w=(rmid+r2)*2+1;

x1=k-round(w/2);
y1=l-round(h/2);
width=w;
height=h;
end