function [t,u] = Trapezi(fun,dfun,t0,T,y0,N)
%Metodo implicito risoluzione ode dei trapezi


%calcolo il passo
h=(T-t0)/N;

t=t0:h:T;

u(1)=y0;
tol=1e-5;
kmax=10;

%dfun nel mio caso è -40

for n=1:N
    %%applico il metodo di Newton
    F=@(z) z-u(n)-h/2*(feval(fun,t(n),u(n))+feval(fun,t(n+1),z)); %%F_EE
    dF=@(z) 1-h/2*feval(dfun,t(n+1),z);
    z0=u(n);
    k=1;d0=1;
    %%metodo di Newton ad ogni step temporale
    while (abs(d0)>tol) & (k<=kmax)%questo per mettere un numero massimo di iterazioni
        d0=-F(z0)/dF(z0);
        z1=z0+d0;
        k=k+1;
        z0=z1;
    end
    u(n+1)=z1;
end
t=t(:); u=u(:);

