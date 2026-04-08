function [t,u] = Trapezisys(fun,Jfun,t0,T,y0,N,epsilon)
%%Metodo di Eulero implicito per sistemi ode
%%y'=f(t,y) f:[t0,T]*R^m->R^m, m>1 sistemi
%%input:
%%fun=M file di tipo function per la f del prob di cauchy
%%Jfun=matrice jacobiana per fun (per il metodo di Newton)


%%output:

h=(T-t0)/N;%passo di discretizzazione

t=t0:h:T;%griglia temporale
tol=1e-5;
kmax=10;
u(:,1)=y0;% la m è la lungh di y0
m=length(y0);I=eye(m);

for n=1:N
    %%applico il metodo di Newton per sistemi ad ogni step
    F=@(z) z-u(:,n)-h/2*(feval(fun,t(n),u(:,n),epsilon)+feval(fun,t(n+1),z,epsilon)); %%F_EI
    JF=@(z) I-h/2*feval(Jfun,t(n+1),z,epsilon);
    z0=u(:,n);
    k=1;d0=1;
    %%metodo di Newton ad ogni step temporale
    while (norm(d0)>tol) & (k<=kmax)%questo per mettere un numero massimo di iterazioni
        d0=-JF(z0)\F(z0);%sistema lineare m*m con metodo diretto LU
        z1=z0+d0;
        k=k+1;
        z0=z1;
    end
    u(:,n+1)=z1;%% valore approssimato da Newton per sistemi
end
t=t(:); 


