# Numerical-Methods-for-ODE
# Metodi Numerici per Equazioni Differenziali Ordinarie (ODE)

[cite_start]Questo repository contiene l'implementazione in MATLAB e l'analisi di convergenza di vari metodi numerici per la risoluzione di Problemi di Cauchy (PdC)[cite: 2185]. [cite_start]Il progetto confronta algoritmi espliciti e impliciti, analizzandone stabilità, ordine di convergenza ed efficienza su problemi stiff e non stiff[cite: 2192, 2721, 2883].

## 🧮 Metodi Implementati (Da zero)
[cite_start]Nella cartella `Solvers` sono presenti le implementazioni *from scratch* dei seguenti metodi[cite: 2629]:
* [cite_start]**Metodi Espliciti:** Eulero Esplicito (EE) [cite: 2655][cite_start], Runge-Kutta 4 (RK4) [cite: 2709][cite_start], Adams-Bashforth 4 (AB4)[cite: 2665].
* [cite_start]**Metodi Impliciti (con Newton):** Trapezi [cite: 2630][cite_start], Adams-Moulton 3 (AM3)[cite: 2681].
* [cite_start]**Sistemi:** RK4sys [cite: 4198][cite_start], Trapezisys[cite: 4197].

[cite_start]Questi metodi vengono confrontati con i solutori integrati di MATLAB (`ode45` e `ode15s`)[cite: 2185].

---

## 🔬 Esercizio 1: Problema Lineare Stiff
[cite_start]Analisi di un problema fortemente stiff caratterizzato da $\lambda = -40$[cite: 2516].
* **Metodi Espliciti:** Mostrano forti limitazioni sul passo $h$. [cite_start]Il metodo RK4, nonostante un'iniziale divergenza dovuta alla sua limitata regione di assoluta stabilità, si dimostra il migliore tra gli espliciti grazie al suo ordine $p=4$, raggiungendo un errore $\epsilon_{max} \approx 4.3 \cdot 10^{-12}$ a $k=15$ iterazioni[cite: 2195]. [cite_start]L'ordine di Eulero Esplicito ($p=1$) risulta invece troppo basso per competere[cite: 2195, 2197].
* [cite_start]**Metodi Impliciti:** Il metodo dei Trapezi si fa notare per la sua A-stabilità, garantendo la convergenza fin dalla prima iterazione senza esplodere[cite: 2190, 2195].
* [cite_start]**Vincitore Assoluto:** `ode15s` si conferma il metodo ideale per problemi stiff, adattando il passo senza restrizioni di stabilità e ottenendo ottimi risultati in sole 76 iterazioni[cite: 2194].

![Confronto Errori RK4](inserisci_qui_il_percorso_immagine_RK4.png)

---

## 📈 Esercizio 2: Problema Non Lineare
[cite_start]Su questo problema, meno stiff del precedente, RK4 performa eccellentemente garantendo un $\epsilon_{max} \approx 1.3 \cdot 10^{-13}$ a $k=15$ iterazioni[cite: 2721].
[cite_start]Tuttavia, `ode15s` si dimostra ancora una volta superiore in efficienza, risolvendo il problema in soli 39 passi con un errore di $\epsilon_{max} \approx 7.2 \cdot 10^{-4}$, mentre RK4 necessita di almeno 64 passi per iniziare a convergere[cite: 2722].

---

## 🌪️ Esercizio 3: Il Sistema di Davis
[cite_start]Questo esercizio analizza un sistema di equazioni differenziali in cui la componente $y_1(t)$ presenta una rigidità dipendente da un parametro $\epsilon$, mentre $y_2(t)$ rimane non stiff (un semplice esponenziale)[cite: 2883].

Sono stati testati tre regimi di stiffness:
1. **Bassa Stiffness ($\epsilon = 1$):** I metodi espliciti dominano. [cite_start]RK4 risulta il migliore, ottenendo un $\epsilon_{max} \approx 2 \cdot 10^{-6}$ in appena 32 passi[cite: 2886].
2. **Stiffness Media ($\epsilon = 10^{-2}$):** Il problema si irrigidisce. [cite_start]`ode15s` diventa il metodo più efficiente (58 passi per un $\epsilon_{max} \approx 1.6 \cdot 10^{-3}$)[cite: 3298]. [cite_start]RK4 e Trapezi richiedono oltre 512 passi per ottenere precisioni accettabili[cite: 3299].
3. [cite_start]**Alta Stiffness ($\epsilon = 10^{-4}$):** La soluzione $y_1(t)$ subisce variazioni repentine, facendo esplodere gli errori iniziali[cite: 3721]. [cite_start]L'errore globale (norma delle due soluzioni) è totalmente dominato da $y_1(t)$[cite: 3725]. [cite_start]L'efficienza di `ode15s` è insuperabile (68 passi), mentre `ode45` arranca richiedendo ben 60289 passi per non perdere la stabilità[cite: 3722, 3723].

![Spazio delle Fasi Sistema di Davis](inserisci_qui_il_percorso_immagine_Davis.png)

---
*Progetto realizzato per il corso di Metodi Numerici.*
