# Numerical-Methods-for-ODE
# Metodi Numerici per Equazioni Differenziali Ordinarie (ODE)

Questo repository contiene l'implementazione in MATLAB e l'analisi di convergenza di vari metodi numerici per la risoluzione di Problemi di Cauchy (PdC). Il progetto confronta algoritmi espliciti e impliciti, analizzandone stabilità, ordine di convergenza ed efficienza su problemi stiff e non stiff.

## 🧮 Metodi Implementati (Da zero)
Nella cartella `Solvers` sono presenti le implementazioni *from scratch* dei seguenti metodi numerici:
* **Metodi Espliciti:** Eulero Esplicito (EE), Runge-Kutta 4 (RK4), Adams-Bashforth 4 (AB4).
* **Metodi Impliciti (con Newton):** Trapezi, Adams-Moulton 3 (AM3).
* **Sistemi:** RK4sys, Trapezisys.

Questi metodi scritti a mano vengono sistematicamente confrontati con i solutori integrati e ottimizzati di MATLAB (`ode45` e `ode15s`).

---

## 🔬 Esercizio 1: Problema Lineare Stiff
Analisi di un problema fortemente stiff caratterizzato da $\lambda = -40$.
* **Metodi Espliciti:** Mostrano forti limitazioni sul passo $h$. Il metodo RK4, nonostante un'iniziale divergenza dovuta alla sua limitata regione di assoluta stabilità, si dimostra il migliore tra gli espliciti grazie al suo ordine teorico $p=4$, raggiungendo un errore $\epsilon_{max} \approx 4.3 \cdot 10^{-12}$ a $k=15$ iterazioni. L'ordine di Eulero Esplicito ($p=1$) risulta invece troppo basso per competere.
* **Metodi Impliciti:** Il metodo dei Trapezi si fa notare per la sua A-stabilità, garantendo la convergenza fin dalla prima iterazione senza esplodere, seppur con un ordine $p=2$.
* **Vincitore Assoluto:** `ode15s` si conferma il metodo ideale per i problemi stiff, adattando il passo dinamicamente senza restrizioni di stabilità e ottenendo ottimi risultati in sole 76 iterazioni.

![Confronto Errori Esercizio 1](inserisci_qui_il_nome_della_foto_1.png)

---

## 📈 Esercizio 2: Problema Non Lineare
Su questo problema, meno stiff del precedente, RK4 performa eccellentemente garantendo un $\epsilon_{max} \approx 1.3 \cdot 10^{-13}$ a $k=15$ iterazioni. 
Risulta addirittura più accurato di `ode45` (che impiega 273 passi con un errore di $\epsilon_{max} \approx 1 \cdot 10^{-3}$). Tuttavia, `ode15s` si dimostra ancora una volta superiore in pura efficienza computazionale, risolvendo il problema in soli 39 passi.

![Soluzione Esercizio 2](inserisci_qui_il_nome_della_foto_2.png)

---

## 🌪️ Esercizio 3: Il Sistema di Davis
Questo esercizio analizza un sistema di equazioni differenziali accoppiate in cui la componente $y_1(t)$ presenta una rigidità fortemente dipendente da un parametro $\epsilon$, mentre $y_2(t)$ rimane non stiff (comportandosi come un semplice esponenziale).

Sono stati testati tre regimi di stiffness:
1. **Bassa Stiffness ($\epsilon = 1$):** I metodi espliciti dominano. RK4 risulta il migliore in assoluto, ottenendo un $\epsilon_{max} \approx 2 \cdot 10^{-6}$ in appena 32 passi.
2. **Stiffness Media ($\epsilon = 10^{-2}$):** Il problema si irrigidisce notevolmente. `ode15s` diventa il metodo più efficiente (58 passi per un $\epsilon_{max} \approx 1.6 \cdot 10^{-3}$). I metodi a passo fisso come RK4 e Trapezi richiedono oltre 512 passi per ottenere precisioni accettabili.
3. **Alta Stiffness ($\epsilon = 10^{-4}$):** La soluzione $y_1(t)$ subisce variazioni repentine (strato limite), facendo esplodere gli errori iniziali se il passo non è sufficientemente piccolo. L'errore globale (calcolato come norma delle due soluzioni) è totalmente dominato da $y_1(t)$. L'efficienza del solutore implicito `ode15s` è insuperabile (impiega solo 68 passi), mentre `ode45` (esplicito) arranca tremendamente, richiedendo ben 60.289 passi per non perdere la stabilità numerica.

![Spazio delle Fasi Sistema di Davis](inserisci_qui_il_nome_della_foto_3.png)

---

## 🔍 Analisi di Convergenza e Ordine dei Metodi
Oltre alla risoluzione dei problemi specifici, il repository include script dedicati allo studio rigoroso della convergenza asintotica dei metodi implementati.

Lo script di test permette di valutare il risolutore su problemi di Cauchy di difficoltà crescente. Il codice effettua un raffinamento spaziale dimezzando iterativamente il passo di integrazione $h$ ($N_k = 2^k$) e calcola l'ordine di convergenza empirico $p$ tramite la formula:

$$p = \log_2 \left( \frac{E(h)}{E(h/2)} \right)$$

I risultati ottenuti confermano empiricamente l'ordine teorico dei metodi trattati (es. $p \approx 1$ per i metodi di Eulero, $p \approx 4$ per Runge-Kutta).

---
*Progetto realizzato per il corso di Metodi Numerici per ODE.*
