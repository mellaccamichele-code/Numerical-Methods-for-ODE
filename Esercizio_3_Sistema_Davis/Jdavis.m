function J=Jdavis(t,y,epsilon)

y1=y(1);
y2=y(2);

f1y1=-1/epsilon;
f1y2=1./(epsilon*(1+y2).^2)-(1-y2)./(1+y2)^3;
f2y1=0;
f2y2=-1;

J=[f1y1,f1y2;f2y1,f2y2];

end