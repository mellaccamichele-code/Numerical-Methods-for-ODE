function [t,u] = EuleroE(fun,t0,T,y0,N)
%%Metodo di Eulero esplicito per ode scalari
%%input:

%%output:

h=(T-t0)/N;
%creo il vettore deu tempi che parte da t0 fino a T con passo h
t=t0:h:T;


u(1)=y0;
for n=1:N
    %feval prende fun e ci mette dentro t e u
    fn=feval(fun,t(n),u(n));
    u(n+1)=u(n)+h*fn;
end
t=t(:); u=u(:);


