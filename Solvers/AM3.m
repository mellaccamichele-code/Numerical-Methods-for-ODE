function [t,u] = AM3(fun,dfun,t0,T,y0,N)
%Metodo di adams-moulton di passo k=3 e p=4 (k+1) implicito per ode scalari
%a k=3 passi, di ordine di converegenza (teorico) p=4


h = (T - t0) / N;
t = t0:h:T; 
tol=1e-5;
kmax=10;

%%parto con RK di ordine p=3, genero u1,u2
[unused,uu]=RK4(fun,t0,t0+2*h,y0,2);
%% uu =[u0,u1,u2]
u(1:3)=uu;%metto i valori trovati in uu in u vettore della soluzione finale
%%
for n=1:N-2 %N-k-1, il ciclo parte da un,un+1,un+2 calcola un+3
    fn=feval(fun,t(n),u(n));
    fn1=feval(fun,t(n+1),u(n+1));
    fn2=feval(fun,t(n+2),u(n+2));

    %applico il metodo di Newton per il u(n+3)
    F=@(z) z-u(n+2)-h/24*(fn-5*fn1+19*fn2+9*feval(fun,t(n+3),z)); %%F_EE
    dF=@(z) 1-h/24*9*feval(dfun,t(n+3),z);%dfun=-40
    z_old=u(n+2);%prendiamo l'ultimo valore noto:u(n+2)
    k=1;d0=1;
    %%metodo di Newton ad ogni step temporale
    while (abs(d0)>tol) & (k<=kmax)%questo per mettere un numero massimo di iterazioni
        d0=-F(z_old)/dF(z_old);
        z_new=z_old+d0;%cioe z_old-F(z_old)/dF(z_old)
        k=k+1;
        z_old=z_new;
    end
    u(n+3)=z_new;

end
t=t(:);u=u(:);
