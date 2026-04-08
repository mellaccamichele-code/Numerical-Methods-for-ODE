%%
% 
%   for x = 1:10
%       disp(x)
%   end
% 
clear all

t0=0;T=5;y0=1; lam=-40;
f=@(t,y) lam*y;
df=@(t,y)-40;
ye=@(t) y0*exp(lam*(t-t0));

figure(1); clf;  

nfig=1;

kmax=15;
for k=1:kmax
    
    N(k)=2^k;

    [tam,uam]=AM3(f,df,t0,T,y0,N(k));
    erram=abs(ye(tam)-uam);
    figure (nfig), subplot(121), plot(tam,uam,'.-'),xlabel('t'),ylabel('y(t)'), title('AM3'),    hold on

    if k>=2
    [tab,uab]=AB4(f,t0,T,y0,N(k));
    errab=abs(ye(tab)-uab);
    figure (nfig+1), subplot(121), plot(tab,uab,'.-'),xlabel('t'),ylabel('y(t)'), title('AB4'),    hold on

    end

    [tee,uee]=EuleroE(f,t0,T,y0,N(k));
    erree=abs(ye(tee)-uee);
    figure (nfig+2), subplot(121), plot(tee,uee,'.-'),xlabel('t'),ylabel('y(t)'), title('EE'),    hold on

    
    [trk,urk]=RK4(f,t0,T,y0,N(k));
    errrk=abs(ye(trk)-urk);
    figure (nfig+3), subplot(121), plot(trk,urk,'.-'),xlabel('t'),ylabel('y(t)'), title('RK4'),    hold on


    [tt,ut]=Trapezi(f,df,t0,T,y0,N(k));
    errt=abs(ye(tt)-ut);
    figure (nfig+4), subplot(121), plot(tt,ut,'.-'),xlabel('t'),ylabel('y(t)'), title('Trapezi'),   hold on




  

    Eam(k)=max(erram);
    if k>=2
    Eab(k)=max(errab);
    end
    Eee(k)=max(erree);
    Erk(k)=max(errrk);
    Et(k)=max(errt);


    figure(nfig), subplot(122), semilogy(tam,erram,'.-'),xlabel('t'),ylabel('err'), title('e_n'), hold on
    if k>=2
    figure(nfig+1), subplot(122), semilogy(tab,errab,'.-'),xlabel('t'),ylabel('err'), title('e_n'), hold on
    end
    figure(nfig+2), subplot(122), semilogy(tee,erree,'.-'),xlabel('t'),ylabel('err'), title('e_n'), hold on
    figure(nfig+3), subplot(122), semilogy(trk,errrk,'.-'),xlabel('t'),xlabel('err'), title('e_n'), hold on
    figure(nfig+4), subplot(122), semilogy(tt,errt,'.-'),xlabel('t'),xlabel('err'), title('e_n'), hold on

  
    if k>1
        pam(k)=log2(Eam(k-1)/Eam(k));
        pab(k)=log2(Eab(k-1)/Eab(k));
        pee(k)=log2(Eee(k-1)/Eee(k));
        prk(k)=log2(Erk(k-1)/Erk(k));
        pt(k)=log2(Et(k-1)/Et(k));
       
    end
    pause(1);
end

[t15s,u15s]=ode15s(f,[t0,T],y0);
 err15s=abs(ye(t15s)-u15s);
    E15s=max(err15s);
   
    figure (nfig+5),subplot(121), plot(t15s,u15s,'.-'),xlabel('t'),ylabel('y(t)'),   title('ode15s'),    hold on
    figure (nfig+5),subplot(122), semilogy(t15s,err15s,'.-'),xlabel('t'),ylabel('err'),   title('e_n'),    hold on

[t45,u45]=ode45(f,[t0,T],y0);
err45=abs(ye(t45)-u45);
    E45=max(err45);
   
    figure (nfig+6),subplot(121), plot(t45,u45,'.-'),xlabel('t'),ylabel('y(t)'),   title('ode45'),    hold on
    figure (nfig+6),subplot(122), semilogy(t45,err45,'.-'),xlabel('t'),ylabel('err'),   title('e_n'),    hold on


    h=(T-t0)./N;

    disp('metodo AM3')
disp('passo h ------ errore max ----- p(ordine)----N' )
ris=[h(:),Eam(:),pam(:),N(:)];
format long g
disp(ris)

disp('metodo AB4')
disp('passo h ------ errore max ----- p(ordine)----N')
ris=[h(:),Eab(:),pab(:),N(:)];
format long g
disp(ris)

disp('metodo EE')
disp('passo h ------ errore max ----- p(ordine)----N')
ris=[h(:),Eee(:),pee(:),N(:)];
format long g
disp(ris)

disp('metodo RK4')
disp('passo h ------ errore max ----- p(ordine)----N')
ris=[h(:),Erk(:),prk(:),N(:)];
format long g
disp(ris)

disp('metodo Trapezi')
disp('passo h ------ errore max ----- p(ordine)----N')
ris=[h(:),Et(:),pt(:),N(:)];
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