# bc_inference

Questa repository contiene due programmi,
uno per generare dati da un modello noto di meccanica statistica (Blume Capel) e l'altro per inferire sulle interazioni utilizzare a partire dai dati generati.
Il modello utilizzato Blume-Capel, che è una generalizzazione del modello di Ising. La sua hamiltoniana è data da:
H(s) = - \sum_{<ij>}J_{ij} s_i s_j +\mu \sum_{i} s_i^2

La somma <ij> è  da intendersi come somma sui nodi di un grafo. Ad esempio in un reticolo 2D i nodi connessi al nodo s0 sono s1, s2, s3, s4:
![3-s2 0-B9780123850300000177-f17-11-9780123850300](https://user-images.githubusercontent.com/86679938/235530157-4250a70b-30cb-4770-80b2-85e77162eedb.jpg)

Il programma fornito è in grado di simulare il modello su un reticolo bidimensionale oppure su un grafo Erdos-Renyi G(n,p) dove p è la probabilità di connessione di ogni coppia di nodi.

Il programma di inferenza prende i dati esportati e tenta di riprodurre la matrice J_{ij} delle interazioni

# Installazione
Per poter utilizzare i due programmi è necessario avere la libreria GSL - GNU Scientific Library.
Ci sono diversi modi per ottenerla. La più semplice, per i sistemi operativi Debian base (e.g. Ubuntu), è quella di digitare sul terminale:

sudo apt-get install libgsl-dev

Altrimenti si può installare manualmente il pacchetto direttamente dal sito: https://mirror.kumi.systems/gnu/gsl/, scaricando il file "gsl-latest.tar.gz".
In questo caso si devono leggere le istruzioni all'interno del file INSTALL presente dentro il pacchetto.

# Come utilizzare
Il programma per generare dati prende come argomenti i valori: L, J, mu, T, p, MCS, count
L è la dimensione del reticolo se si esegue il programma su reticolo, altrimenti N=L*L rappresenta il numero di spins.
J è la costante di interazione (di solito lasciata = 1, ma uno può provare con diversi valori).
mu è il potenziale chimico che favorisce la presenza si spin nulli.
p è la probabilità per impostare il grafo G(n,p)
MCS è il numero di steps monte carlo 
count è una variabile per numerare i file esportati

Oltre a questi parametri si possono modificare due #define all'interno del programma:
TYPE, 1 indica il reticolo, 0 il grafo
THERM, MCS iniziali per termalizzare in modo ragionevole (default 100000 dovrebbero essere sufficienti)

Il programma di inferenza invece riceve: L, MCS, Tmin, Tmax, dT, meas
Tmin, Tmax, dT sono rispettivamente la temperatura minima, massima, passo di incremento
meas è il numero di misure fatte

I file vengono generati dentro la directory /outputs.
