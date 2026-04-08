function [t,u] = AB4(fun,t0,T,y0,N)
%Metodo di adams-bashforth esplicito per ode scalari
%a k=4 passi, di ordine di converegenza (teorico) p=4


h = (T - t0) / N;
t = t0:h:T; 




%%parto con RK di ordine p=4, genero u1,u2,u3
[unused,uu]=RK4(fun,t0,t0+3*h,y0,3);
%% uu =[u0,u1,u2,u3]
u(1:4)=uu;%metto i valori trovati in uu in u vettore della soluzione finale
%%
for n=1:N-3 %N-k-1
    fn=feval(fun,t(n),u(n));
    fn1=feval(fun,t(n+1),u(n+1));
    fn2=feval(fun,t(n+2),u(n+2));
    fn3=feval(fun,t(n+3),u(n+3));
    u(n+4)=u(n+3)+(h/24)*(-9*fn+37*fn1-59*fn2+55*fn3);
end
t=t(:);u=u(:);
