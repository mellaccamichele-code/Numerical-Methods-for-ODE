% script per studiare la convergenza di RKE di ordine 4
clear all
ode=input('Eulero E esempio di PdC ode=')
switch(ode)
    %% esempi con soluzione esatta magari di difficoltà crescente
    case 1 % problema lineare
        t0=0;T=1;y0=1;
        f=@(t,y) -y+t;
        df=@(t,y) -1;
        ye=@(t) 2*exp(-t)+t-1;
    case 2 % problema NON lineare
        t0=0;T=2;y0=0;
        f=@(t,y) cos(2*y);
        df=@(t,y)-2*sin(2*y);
        ye=@(t)0.5*asin((exp(4*t)-1)./(exp(4*t)+1));
    case 3 %prb test
        t0=0;T=1;y0=1; lam=-100;%-100 %-10
        f=@(t,y) lam*y;
        df=@(t,y)-100;
        ye=@(t) y0*exp(lam*(t-t0));
end
nfig=100*ode;
kmax=10;
for k=1:kmax
    N(k)=2^k;
    [t,u]=AM3(f,df,t0,T,y0,N(k),nfig);
    figure(nfig),hold on
    err=abs(ye(t)-u);
    E(k)=max(err);
    figure(nfig+1), semilogy(t,err,'.-'), hold on, title('errore assoluto e_n')
    pause(1)
    %% stimo la velocità di convergenza
    if k>1
        p(k)=log2(E(k-1)/E(k));
    end
end
tt=linspace(t0,T,500);
figure(nfig), plot(tt,ye(tt),'r.-','LineWidth',2), xlabel('t'),%linea più spessa
hold off
figure(nfig+1), hold off
%% tabella dei risultati 
h=(T-t0)./N;
disp('metodo RK4')
disp('passo h ------ errore max ----- p(ordine)----' )
ris=[h(:),E(:),p(:)];
format long g
disp(ris)