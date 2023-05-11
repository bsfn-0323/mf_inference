#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

class Ising2D {
    private:
        int L; //Lato del reticolo
        double J; //Costante di accoppiamento
        double H; //Campo magnetico esterno
        double T; //Temperatura
        double beta; //Inverso di T
        int N; //Numero totale di spin

        double magn; //Magnetizzazione intensiva
        double energy; //Energia

        mt19937 rng; // Generatore di numeri casuali
        uniform_real_distribution<double> randreal; // Generatore di numeri casuali reali tra 0 e 1

        int *spins; //Array degli spin

    public:
        Ising2D(int _L, double _J, double _H, double _T, bool hot_start = true) : L(_L), J(_J), H(_H), T(_T), beta(1.0/T), N(L*L), randreal(0.0, 1.0), rng(random_device{}()) {
            spins = new int[N];
            //Riempiamo l'array degli spin randomicamente
            if (hot_start){
                for(int i=0; i<N; i++){
                    spins[i] = floor(randreal(rng) * 2) * 2 - 1;  
                }

            //Riempiamo l'array degli spin assegnando tutti +1
            } else {
                for(int i=0; i<N; i++){
                    spins[i] = 1; 
                }
            }
            calculate_magnetization();
            calculate_energy();
        }

        ~Ising2D(){
            delete[] spins;
        }

        void montecarlo_sweep(){
            int selected_spin;
            double deltaE;
            for (int k=0; k<N; k++){
                selected_spin = floor(randreal(rng) * N);
                deltaE = 2*J*spins[selected_spin]*
                    (spins[(selected_spin+1)%N] + spins[(selected_spin-1+N)%N] + 
                    spins[(selected_spin+L)%N] + spins[(selected_spin-L+N)%N]) + 
                    2*H*spins[selected_spin];
                if(deltaE<0 || randreal(rng) < exp(-beta*deltaE)){
                    spins[selected_spin] *= -1;
                }
            }
        }

        double calculate_magnetization(){
            magn = 0.0;
            for(int i=0; i<N; i++){
                magn += spins[i];
            }
            magn /= N;
            return magn;
        }

        double calculate_energy(){
            energy = 0.0;
            for(int i=0; i<N; i++){
                energy -= J * spins[i] * (spins[(i+1)%N] + spins[(i+L)%N]); //Termine di accoppiamento
                energy -= H * spins[i]; //Termine di campo magnetico esterno
            }
            energy /= N;
            return energy;
        }

        double get_energy(){
            return energy;
        }

        double get_magn(){
            return magn;
        }

        int* get_spins(){
            return spins;
        }
        

};

void progressBar(int progress, int total) {
    string bar;
    float percentage = (float)progress / (float)total;
    int barWidth = 70;

    for (int i = 0; i < barWidth * percentage; i++) {
        bar += "#";
    }

    cout << "[" << bar;
    for (int i = 0; i < barWidth - bar.length(); i++) {
        cout << " ";
    }
    cout << "] " << int(percentage * 100.0) << " %\r";
    cout.flush();
}

int main(){
    const int NUM_TEMPS = 1;
    const int NUM_SWEEPS = 1000;
    const int L = 50;
    const double J = 1;
    const double H = 0;

    double temperatures[NUM_TEMPS];
    double magnetizations[NUM_TEMPS];
    double energies[NUM_TEMPS];

    // Generazione delle temperature
    for(int i=0; i<NUM_TEMPS; i++){
        temperatures[i] = i*0.4 + 0.1;
    }

    // Ciclo sulle temperature
    progressBar(0, NUM_TEMPS);
    for(int t=0; t<NUM_TEMPS; t++){
        double T = temperatures[t];
        Ising2D lattice(L, J, H, T, false); // Inizializzazione del reticolo random
        for(int i=0; i<NUM_SWEEPS; i++){
            lattice.montecarlo_sweep(); // Sweep di Montecarlo
        }
        lattice.calculate_magnetization();
        lattice.calculate_energy();
        magnetizations[t] = lattice.get_magn(); // Salvataggio magnetizzazione finale
        energies[t] = lattice.get_energy(); // Salvataggio energia finale
        progressBar(t+1, NUM_TEMPS);
    }

    // Stampa a schermo delle magnetizzazioni ed energie
    for(int i=0; i<NUM_TEMPS; i++){
        cout << "T = " << temperatures[i] << " | Magnetization = " << magnetizations[i] << " | Energy = " << energies[i] << endl;
    }
    
    // Salvataggio dei dati in un file di testo
    ofstream outfile("Ising2D.txt");
    if (outfile.is_open()){
        for(int i=0; i<NUM_TEMPS; i++){
            outfile << temperatures[i] << " " << magnetizations[i] << " " << energies[i] << endl;
        }
        outfile.close();
    } else {
        cout << "Impossibile aprire il file di output." << endl;
    }
}