#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

class Ising2D {
    private:
        int L; //Lato del reticolo
        double** J; //Costante di accoppiamento
        double* H; //Campo magnetico esterno
        double T; //Temperatura
        double beta; //Inverso di T
        int N; //Numero totale di spin

        double magn; //Magnetizzazione intensiva
        double energy; //Energia

        mt19937 rng; // Generatore di numeri casuali
        uniform_real_distribution<double> randreal; // Generatore di numeri casuali reali tra 0 e 1

        int *spins; //Array degli spin

        void create_couplings(){
            J = new double*[N];
            for (int k=0; k<N; k++){
                J[k] = new double[4];
            }
            for(int k=0; k<N; k++){
                J[k][0] = randreal(rng) * 2 - 1; //0: Spin a destra (k+1)
                J[k][2] = randreal(rng) * 2 - 1; //2: Spin in basso (k+L)

                J[(k+1)%N][1] = J[k][0]; //1: Spin a sinistra (k-1)
                J[(k+L)%N][3] = J[k][2]; //3: Spin in alto (k-L)
            }
        }

        void create_fields(){
            H = new double[N];
            for (int k=0; k<N; k++){
                H[k] = 0.;
            }
        }

    public:
        Ising2D(int _L, double _T, bool hot_start = true) : L(_L), T(_T), beta(1.0/T), N(L*L), randreal(0.0, 1.0), rng(random_device{}()) {
            spins = new int[N];
            create_couplings();
            create_fields();
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
            for(int k=0; k<N; k++){
                delete[] J[k];
            }
            delete[] J;
            delete[] H;
        }

        void montecarlo_sweep(){
            int selected_spin;
            double deltaE;
            for (int k=0; k<N; k++){
                selected_spin = floor(randreal(rng) * N);
                deltaE = 2*spins[selected_spin]*
                    (spins[(selected_spin+1)%N]*J[selected_spin][0] + spins[(selected_spin-1+N)%N]*J[selected_spin][1] + 
                    spins[(selected_spin+L)%N]*J[selected_spin][2] + spins[(selected_spin-L+N)%N]*J[selected_spin][3]) + 
                    2*H[selected_spin]*spins[selected_spin];
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
                energy -= spins[i] * (spins[(i+1)%N]*J[i][0] + spins[(i+L)%N]*J[i][2]); //Termine di accoppiamento
                energy -= H[i] * spins[i]; //Termine di campo magnetico esterno
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
    const int NUM_TEMPS = 20;
    const int NUM_SWEEPS = 10;
    const int L = 100;

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
        Ising2D lattice(L, T, true); // Inizializzazione del reticolo random
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
    ofstream outfile("Ising2D_dis.txt");
    if (outfile.is_open()){
        for(int i=0; i<NUM_TEMPS; i++){
            outfile << temperatures[i] << " " << magnetizations[i] << " " << energies[i] << endl;
        }
        outfile.close();
    } else {
        cout << "Impossibile aprire il file di output." << endl;
    }
}