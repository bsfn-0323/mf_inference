using namespace std;

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>

using namespace std;

class IsingER {
    private:
        int N; //Numero totale di spin
        double p; //Probabilità di formazione degli archi
        vector<vector<int>> graph; //Liste di adiacenza
        vector<vector<double>> J; //Costante di accoppiamento
        vector<double> H; //Campo magnetico esterno
        double T; //Temperatura
        double beta; //Inverso di T

        double magn; //Magnetizzazione intensiva
        double energy; //Energia

        mt19937 rng; // Generatore di numeri casuali
        uniform_real_distribution<double> randreal; // Generatore di numeri casuali reali tra 0 e 1

        vector<int> spins; //Array degli spin

        void create_graph(){ //Creazione del grafo random
            graph.resize(N);
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    if (randreal(rng) < p) {
                        graph[i].push_back(j);
                        graph[j].push_back(i);
                    }
                }
            }
        }

        void create_couplings(){ //Creazione dei couplings
            J.resize(N, vector<double>(N, 0.0));
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < graph[i].size(); j++) {
                    J[i][graph[i][j]] = 1.0;
                }
            }
        }

        void create_fields(){
            H.resize(N, 0.0);
        }

    public:
        IsingER(int _N, double _p, double _T, bool hot_start = true) : N(_N), p(_p), T(_T), beta(1.0/T), randreal(0.0, 1.0), rng(random_device{}()) {
            spins.resize(N, 0);
            create_graph();
            create_couplings();
            create_fields();
            //Riempiamo l'array degli spin randomicamente
            if (hot_start) {
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

        void montecarlo_sweep(){
            int selected_spin;
            int interacting_spin;
            double deltaE;
            for (int k=0; k<N; k++){
                selected_spin = floor(randreal(rng) * N);
                deltaE = 0.;
                for(int j=0; j<graph[selected_spin].size(); j++){
                    interacting_spin = graph[selected_spin][j];
                    deltaE += 2*spins[selected_spin]*spins[interacting_spin]*J[selected_spin][interacting_spin];
                }
                deltaE += 2*H[selected_spin]*spins[selected_spin];
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
                for(int j=i+1; j<N; j++){
                    energy -= spins[i]*spins[j]*J[i][j];
                }
                energy -= spins[i]*H[i];
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

        vector<int> get_spins(){
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
    const int NUM_TEMPS = 20;
    const int NUM_SWEEPS = 10000;
    const int N = 100;
    const double p = 0.001;

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
        IsingER lattice(N, p, T, true); // Inizializzazione del reticolo random
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
    ofstream outfile("IsingER.txt");
    if (outfile.is_open()){
        for(int i=0; i<NUM_TEMPS; i++){
            outfile << temperatures[i] << " " << magnetizations[i] << " " << energies[i] << endl;
        }
        outfile.close();
    } else {
        cout << "Impossibile aprire il file di output." << endl;
    }
}