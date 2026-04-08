function dydt=dydt_Davis(t,y,epsilon)
%devo restituire dy/dt il sistema di ode

%input:
%t 
%y, cioè y1=y(1) e y2=y(2)

y1=y(1);
y2=y(2);

f1=(1/epsilon)*(-y1+y2./(1+y2))-y2./(1+y2).^2;
f2=-y2;

dydt=[f1,f2]';% vettore colonna se no mi sfancula

end