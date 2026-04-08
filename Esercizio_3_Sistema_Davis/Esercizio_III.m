clear all


alpha=1.5;
beta=1;
epsilon=[1,10^(-2),10^(-4)];
y0=[alpha,beta];
t0=0;
T=5;

k0=1;
kmax=15;
j=1;
z=13;
tt=linspace(t0,T,1000);

for i=1:3
    
ye1 = @(t) (alpha - beta/(1+beta))*exp(-t/epsilon(i)) + beta./(exp(t)+ beta);
ye2=@(t) beta*exp(-t);


odefun = @(t, y) dydt_Davis(t, y, epsilon(i));%funzione da passare a ode45 e ode15s
[t45,u45]=ode45(odefun,[t0,T],y0);

err45_1=abs(ye1(t45)-u45(:,1));
err45_2=abs(ye2(t45)-u45(:,2));
E45_1(i)=max(err45_1);
E45_2(i)=max(err45_2);

l45(i)=length(t45);

figure(j), subplot(221), plot(t45,u45(:,1),'b-',t45,u45(:,2),'r-',tt,ye1(tt),'k--',tt,ye2(tt),'k--'),xlabel('t'),ylabel('y(t)'), title('ode45');
figure(j), subplot(222), semilogy(t45,err45_1,'b.-',t45,err45_2,'r.-'),xlabel('t'),ylabel('err'), title('e_n');
figure(j), subplot(223), plot(t45(1:end-1), t45(2:end)-t45(1:end-1),'o-');%%creo vettore dei passi sottraendo vettore passi da 1 in poi a vettore dei passi da 2 in poi
figure(z), plot(u45(:,1),u45(:,2),'o-'),xlabel('y_1(t)'),ylabel('y_2(t)'), title ('ode45 spazio delle fasi')%%spazio delle fasi



[t15s,u15s]=ode15s(odefun,[t0,T],y0);
err15s_1=abs(ye1(t15s)-u15s(:,1));
err15s_2=abs(ye2(t15s)-u15s(:,2));
E15s_1(i)=max(err15s_1);
E15s_2(i)=max(err15s_2);

E45_norm(i)  = sqrt(E45_1(i)^2  + E45_2(i)^2);
E15s_norm(i) = sqrt(E15s_1(i)^2 + E15s_2(i)^2);


l15s(i)=length(t15s);

figure(j+1), subplot(221), plot(t15s,u15s(:,1),'b-',t15s,u15s(:,2),'r-',tt,ye1(tt),'k--',tt,ye2(tt),'k--'),xlabel('t'),ylabel('y(t)'), title('ode15s');
figure(j+1), subplot(222), semilogy(t15s,err15s_1,'b.-',t15s,err15s_2,'r.-'),xlabel('t'),ylabel('err'), title('e_n');
figure(j+1), subplot(223), plot(t15s(1:end-1), t15s(2:end)-t15s(1:end-1),'o-');
figure(z+1), plot(u15s(:,1),u15s(:,2),'o-'),xlabel('y_1(t)'),ylabel('y_2(t)'), title ('ode15s spazio delle fasi')%%spazio delle fasi





for k=k0:kmax
N(k)=2^k;
[ttr,ut]=Trapezisys(@dydt_Davis,@Jdavis,t0,T,y0,N(k),epsilon(i));
[trk4,urk4]=RK4sys(@dydt_Davis,t0,T,y0,N(k),epsilon(i));

errt_1=abs(ye1(ttr)-ut(1,:)');
errt_2=abs(ye2(ttr)-ut(2,:)');
    Et_1(i,k)=max(errt_1);1
    Et_2(i,k)=max(errt_2);
    errt_norm = sqrt(errt_1.^2 + errt_2.^2);
    Et_norm(i,k) = max(errt_norm);
    

errrk4_1=abs(ye1(ttr)-urk4(1,:).');
errrk4_2=abs(ye2(ttr)-urk4(2,:).');
    Erk4_1(i,k)=max(errrk4_1);
    Erk4_2(i,k)=max(errrk4_2);
    errrk4_norm = sqrt(errrk4_1.^2 + errrk4_2.^2);
    Erk4_norm(i,k) = max(errrk4_norm);

    
    if k>1
        pt_1(i,k)=log2(Et_1(i,k-1)/Et_1(i,k));
        pt_2(i,k)=log2(Et_2(i,k-1)/Et_2(i,k));
        prk4_1(i,k)=log2(Erk4_1(i,k-1)/Erk4_1(i,k));
        prk4_2(i,k)=log2(Erk4_2(i,k-1)/Erk4_2(i,k));
        pt_norm(i,k)   = log2(Et_norm(i,k-1)/Et_norm(i,k));
        prk4_norm(i,k) = log2(Erk4_norm(i,k-1)/Erk4_norm(i,k));
    end

figure(j+2),subplot(121), plot(ttr,ut(1,:),'b.-',ttr,ut(2,:),'r.-',tt,ye1(tt),'k--',tt,ye2(tt),'k--'),xlabel('t'),ylabel('y(t)'), title('Trapezi'),hold on
figure(j+2),subplot(122), semilogy(ttr,errt_1,'b.-',ttr,errt_2,'r.-'),xlabel('t'),ylabel('err'), title('e_n'),hold on


figure(j+3),subplot(121), plot(trk4,urk4(1,:),'b.-',trk4,urk4(2,:),'r.-',tt,ye1(tt),'k--',tt,ye2(tt),'k--'),xlabel('t'),ylabel('y(t)'), title('RK4'), hold on
figure(j+3),subplot(122), semilogy(trk4,errrk4_1,'b.-',trk4,errrk4_2,'r.-'),xlabel('t'),ylabel('err'), title('e_n'),hold on


if k==kmax
    figure(z+2), plot(ut(1,:),ut(2,:),'o-'),xlabel('y_1(t)'),ylabel('y_2(t)'), title ('Trapezi spazio delle fasi')%%spazio delle fasi
    figure(z+3), plot(urk4(1,:),urk4(2,:),'o-'),xlabel('y_1(t)'),ylabel('y_2(t)'), title ('RK4 spazio delle fasi')%%spazio delle fasi
end

end
j=j+4;
z=z+4;
end

h=(T-t0)./N;

%% tabella dei risultati 
disp('metodo Trapezi con epsilon = 1')
disp('passo h ------ errore max soluzione 1 ------ errore max soluzione 2-----errore in norma ----- p(ordine) per soluzione 1----- p(ordine) per soluzione 2-----p(ordine in norma)----N' )
ris=[h(k0:kmax)',Et_1(1,k0:kmax)',Et_2(1,k0:kmax)',Et_norm(1,k0:kmax)',pt_1(1,k0:kmax)',pt_2(1,k0:kmax)',pt_norm(1,k0:kmax)',N(k0:kmax)'];
format long g
disp(ris)

disp('metodo Trapezi con epsilon = 10^(-2)')
disp('passo h ------ errore max soluzione 1 ------ errore max soluzione 2-----errore in norma ----- p(ordine) per soluzione 1----- p(ordine) per soluzione 2-----p(ordine in norma)----N'  )
ris=[h(k0:kmax)',Et_1(2,k0:kmax)',Et_2(2,k0:kmax)',Et_norm(2,k0:kmax)',pt_1(2,k0:kmax)',pt_2(2,k0:kmax)',pt_norm(2,k0:kmax)',N(k0:kmax)'];
format long g
disp(ris)

disp('metodo Trapezi con epsilon = 10^(-4)')
disp('passo h ------ errore max soluzione 1 ------ errore max soluzione 2-----errore in norma ----- p(ordine) per soluzione 1----- p(ordine) per soluzione 2-----p(ordine in norma)----N'  )
ris=[h(k0:kmax)',Et_1(3,k0:kmax)',Et_2(3,k0:kmax)',Et_norm(3,k0:kmax)',pt_1(3,k0:kmax)',pt_2(3,k0:kmax)',pt_norm(3,k0:kmax)',N(k0:kmax)'];
format long g
disp(ris)

 disp('metodo RK4 con epsilon = 1')
disp('passo h ------ errore max soluzione 1 ------ errore max soluzione 2-----errore in norma ----- p(ordine) per soluzione 1----- p(ordine) per soluzione 2-----p(ordine in norma)----N'  )
ris=[h(k0:kmax)',Erk4_1(1,k0:kmax)',Erk4_2(1,k0:kmax)',Erk4_norm(1,k0:kmax)',prk4_1(1,k0:kmax)',prk4_2(1,k0:kmax)',prk4_norm(1,k0:kmax)',N(k0:kmax)'];
format long g
disp(ris)

disp('metodo RK4 con epsilon = 10^(-2)')
disp('passo h ------ errore max soluzione 1 ------ errore max soluzione 2-----errore in norma ----- p(ordine) per soluzione 1----- p(ordine) per soluzione 2-----p(ordine in norma)----N'  )
ris=[h(k0:kmax)',Erk4_1(2,k0:kmax)',Erk4_2(2,k0:kmax)',Erk4_norm(2,k0:kmax)',prk4_1(2,k0:kmax)',prk4_2(2,k0:kmax)',prk4_norm(2,k0:kmax)',N(k0:kmax)'];
format long g
disp(ris)

disp('metodo RK4 con epsilon = 10^(-4)')
disp('passo h ------ errore max soluzione 1 ------ errore max soluzione 2-----errore in norma ----- p(ordine) per soluzione 1----- p(ordine) per soluzione 2-----p(ordine in norma)----N' )
ris=[h(k0:kmax)',Erk4_1(3,k0:kmax)',Erk4_2(3,k0:kmax)',Erk4_norm(3,k0:kmax)',prk4_1(3,k0:kmax)',prk4_2(3,k0:kmax)',prk4_norm(2,k0:kmax)',N(k0:kmax)'];
format long g
disp(ris)

 


disp('metodo ode15s con epsilon = 1')
disp('errore max soluzione 1 ------ errore max soluzione 2 -----errore in norma------ N' )
ris=[E15s_1(1),E15s_2(1),E15s_norm(1),l15s(1)];
format long g
disp(ris)
% 

disp('metodo ode15s con epsilon = 10^(-2)')
disp('errore max soluzione 1 ------ errore max soluzione 2 -----errore in norma------ N' )
ris=[E15s_1(2),E15s_2(2),E15s_norm(2),l15s(2)];
format long g
disp(ris)
% 

disp('metodo ode15s con epsilon = 10^(-4)')
disp('errore max soluzione 1 ------ errore max soluzione 2 -----errore in norma------ N')
ris=[E15s_1(3),E15s_2(3),E15s_norm(3),l15s(3)];
format long g
disp(ris)

disp('metodo ode45 con epsilon = 1')
disp('errore max soluzione 1 ------ errore max soluzione 2 -----errore in norma------ N' )
ris=[E45_1(1),E45_2(1),E45_norm(1),l45(1)];
format long g
disp(ris)
% 

disp('metodo ode45 con epsilon = 10^(-2)')
disp('errore max soluzione 1 ------ errore max soluzione 2 -----errore in norma------ N' )
ris=[E45_1(2),E45_2(2),E45_norm(2),l45(2)];
format long g
disp(ris)
% 

disp('metodo ode45 con epsilon = 10^(-4)')
disp('errore max soluzione 1 ------ errore max soluzione 2 -----errore in norma------ N' )
ris=[E45_1(3),E45_2(3),E45_norm(3),l45(3)];
format long g
disp(ris)



clear all