function [t,u]=RK4(fun,t0,T,y0,N,epsilon)
%% METODO DI RUNGE-KUTTA 4 PER ODE SCALARI
%% INPUT:

%% OUTPUT:

h=(T-t0)/N;
t=t0:h:T;
u(:,1)=y0;
for n=1:N
    fn=feval(fun,t(n),u(:,n),epsilon);
    %% calcolo i 4 stages 
    K1 = fn;
    K2 = feval(fun, t(n) + h/2, u(:,n) + h * K1/2,epsilon);
    K3 = feval(fun, t(n) + h/2, u(:,n) + h * K2/2,epsilon);
    K4 = feval(fun, t(n) + h, u(:,n) + h * K3,epsilon);
    
    u(:,n+1) = u(:,n) + (h/6)*(K1+2*K2+2*K3+K4);
end
t=t(:);
