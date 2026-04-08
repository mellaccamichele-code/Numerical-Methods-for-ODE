close all

t0=1;T=20;y1=1; 
f=@(t,y) -5.*t.*(y.^2)+(5.*t-1)./(t.^2);
df=@(t,y)-10*t*y;
ye=@(t) 1./t;
nfig=1;
kmax=15;

tt=linspace(t0,T,500);
figure(nfig),subplot(121), plot(tt,ye(tt),'r.-','LineWidth',2), hold on%linea più spessa

for k=6:kmax
    N(k)=2^k;
    [trk,urk]=RK4(f,t0,T,y1,N(k));
    figure(nfig), subplot(121), plot(trk,urk,'b-'),xlabel('t'),ylabel('y(t)'), hold on, title('RK4')

   
  
    
    err=abs(ye(trk)-urk);
    E(k)=max(err);
    figure(nfig), subplot(122), semilogy(trk,err,'.-'), hold on, title('errore assoluto e_n')
    pause(1)
    %% stimo la velocità di convergenza
    if k>1
        p(k)=log2(E(k-1)/E(k));
    end
end

[t15s,u15s]=ode15s(f,[t0,T],y1);
   
    figure (nfig+1), subplot(121), plot(t15s,u15s,'o-'),xlabel('t'),ylabel('y(t)'),   title('ode15s'),    hold on


    err15s=abs(ye(t15s)-u15s);
    E15s=max(err15s);
    figure(nfig+1), subplot(122), semilogy(t15s,err15s,'.-'), hold on, title('e_n')
    pause(1)
    %% stimo la velocità di convergenza

    figure(nfig+1),subplot(121), plot(tt,ye(tt),'r.-','LineWidth',1), hold on%linea più spessa




[t45,u45]=ode45(f,[t0,T],y1);
   
    figure (nfig+2), plot(t45,u45,'o-'),xlabel('t'),ylabel('y(t)'),   title('ode45'),    hold on
        


    err45=abs(ye(t45)-u45);
    E45=max(err45);
    figure(nfig+2), subplot(122), semilogy(t45,err45,'.-'), hold on, title('e_n')
    pause(1)
    %% stimo la velocità di convergenza

    figure(nfig+2),subplot(121), plot(tt,ye(tt),'r.-','LineWidth',1), hold on%linea più spessa

%% tabella dei risultati 
h=(T-t0)./N;
disp('metodo RK4')
disp('passo h ------ errore max ----- p(ordine)----N' )
ris=[h(:),E(:),p(:),N(:)];
format long g
disp(ris)

disp('metodo ode15s')
disp('errore max ----- N' )
ris=[E15s,length(t15s)];
format long g
disp(ris)

disp('metodo ode45')
disp('errore max ----- N' )
ris=[E45,length(t45)];
format long g
disp(ris)

clear all