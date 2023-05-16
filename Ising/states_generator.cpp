#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

//Ising on a 2D lattice
class Ising2D {
    private:
        int L; //Lato del reticolo
        double** J; //Costante di accoppiamento
        double* H; //Campo magnetico esterno
        double T; //Temperatura
        double beta; //Inverso di T
        int N; //Numero totale di spin
        bool is_random; //true se è uno spin glass, false se è ordinato

        double magn; //Magnetizzazione intensiva
        double energy; //Energia

        mt19937 rng; // Generatore di numeri casuali
        uniform_real_distribution<double> randreal; // Generatore di numeri casuali reali tra 0 e 1;
        normal_distribution<double> randgauss;
        vector<int> spins; //Array degli spin

        void create_couplings(bool is_random){
            J = new double*[N];
            for (int k=0; k<N; k++){
                J[k] = new double[4];
            }

            if(is_random){
                for(int k=0; k<N; k++){
                    J[k][0] = randgauss(rng); //0: Spin a destra (k+1)
                    J[k][2] = randgauss(rng); //2: Spin in basso (k+L)
                    J[(k+1)%N][1] = J[k][0]; //1: Spin a sinistra (k-1)
                    J[(k+L)%N][3] = J[k][2]; //3: Spin in alto (k-L)
                }
            } else {
                for(int k=0; k<N; k++){
                    J[k][0] = 1.0; //0: Spin a destra (k+1)
                    J[k][2] = 1.0; //2: Spin in basso (k+L)
                    J[(k+1)%N][1] = J[k][0]; //1: Spin a sinistra (k-1)
                    J[(k+L)%N][3] = J[k][2]; //3: Spin in alto (k-L)
                }

            }
        }

        void create_fields(){
            H = new double[N];
            for (int k=0; k<N; k++){
                H[k] = 0.;
            }
        }

    public:
        Ising2D(int _L, double _T, bool _is_random, bool hot_start = true) : L(_L), T(_T), is_random(_is_random), beta(1.0/T), N(L*L), randreal(0.0, 1.0), rng(random_device{}()), randgauss(0.0, 1.0) {
            spins.resize(N, 0);
            create_couplings(is_random);
            create_fields();
            spin_init(hot_start);
            calculate_magnetization();
            calculate_energy();
        }

        ~Ising2D(){
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

        vector<int> get_spins(){
            return spins;
        }

        vector<vector<double>> get_couplings(){
            vector<vector<double>> couplings(N, vector<double>(N, 0.));
            for(int k=0; k<N; k++){
                couplings[k][(k+1)%N] = J[k][0];
                couplings[k][(k+L)%N] = J[k][2];
                couplings[k][(k-1+N)%N] = J[k][1];
                couplings[k][(k-L+N)%N] = J[k][3];
            }
            return couplings;
        }

        void spin_init(bool hot_start){
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

        }

};

//Ising on a random graph
class IsingER {
    private:
        int N; //Numero totale di spin
        double p; //Probabilità di formazione degli archi
        vector<vector<int>> graph; //Liste di adiacenza
        vector<vector<double>> J; //Costante di accoppiamento
        vector<double> H; //Campo magnetico esterno
        double T; //Temperatura
        double beta; //Inverso di T
        bool is_random; //true se è uno spin glass, false se è ordinato

        double magn; //Magnetizzazione intensiva
        double energy; //Energia

        mt19937 rng; // Generatore di numeri casuali
        uniform_real_distribution<double> randreal; // Generatore di numeri casuali reali tra 0 e 1
        normal_distribution<double> randgauss;
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

        void create_couplings(bool is_random){ //Creazione dei couplings
            J.resize(N, vector<double>(N, 0.0));
            if(is_random){
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < graph[i].size(); j++) {
                        J[i][graph[i][j]] = randgauss(rng);
                        J[graph[i][j]][i] = J[i][graph[i][j]];
                    }
                }
            } else{
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < graph[i].size(); j++) {
                        J[i][graph[i][j]] = 1.0;
                        J[graph[i][j]][i] = J[i][graph[i][j]];
                    }
                }
            }
        }

        void create_fields(){
            H.resize(N, 0.0);
        }

    public:
        IsingER(int _N, double _p, double _T, bool _is_random, bool hot_start = true) : N(_N), p(_p), T(_T), is_random(_is_random), beta(1.0/T), randreal(0.0, 1.0), rng(random_device{}()), randgauss(0.0, 1.0) {
            spins.resize(N, 0);
            create_graph();
            create_couplings(is_random);
            create_fields();
            spin_init(hot_start);
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

        vector<vector<double>> get_couplings(){
            return J;
        }

        void spin_init(bool hot_start){
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

        }
        
};

int main(int argc, char *argv[]){

    //Possiamo scegliere tra Lattice e random graph attivando il corrispondente define

    #define LATTICE_2D
    //#define RANDOM_GRAPH

    #ifdef LATTICE_2D

    if(argc != 4){
        cout << "usage: " << argv[0] << " <T> <L> <is_dis>"<< endl;
        return 1;
    }

    const int NUM_SWEEPS = 200000; //Numero di sweep montecarlo
    const double T = atof(argv[1]); //Temperatura 
    const int L = atoi(argv[2]); //Lato del reticolo
    const bool is_dis = atoi(argv[3]); //true se è uno spin glass
    int N = L * L; //Numero di spin
    
    Ising2D ising(L, T, is_dis);
    vector<vector<double>> couplings = ising.get_couplings();
    stringstream interaction_stream;
    interaction_stream << "./data_lattice/interaction_L" << L << "_T" << fixed << setprecision(2) << T << (is_dis ? "_dis" : "_ord") << ".dat";
    string interaction_filename = interaction_stream.str();
    ofstream outfile1(interaction_filename);

    stringstream config_stream;
    config_stream << "./data_lattice/config_L" << L << "_T" << fixed << setprecision(2) << T << (is_dis ? "_dis" : "_ord") << ".dat";
    string config_filename = config_stream.str();
    ofstream outfile2(config_filename);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            outfile1 << couplings[i][j] << "\t";
        }
        outfile1 << "\n"; 
    }

    for(int j=0; j<NUM_SWEEPS; j++){
        ising.montecarlo_sweep(); // Sweep di Montecarlo
        vector<int> spins = ising.get_spins();
        for (int k = 0; k < N; k++) {
            outfile2 << spins[k] << "\t"; 
        }
        outfile2 << endl;
    }
    return 129;
    #endif

    #ifdef RANDOM_GRAPH

    if(argc != 4){
        cout << "usage: " << argv[0] << " <T> <N> <is_dis>"<< endl;
        return 1;
    }

    const int NUM_SWEEPS = 200000; //Numero di sweep montacarlo
    const double T = atof(argv[1]); //Temperatura
    const int N = atoi(argv[2]); //Numero di spin
    const bool is_dis = atoi(argv[3]); //true se è uno spin glass
    const double p = 0.04; //Probabilità di formazione di un legame

    IsingER ising(N, p, T, is_dis);
    vector<vector<double>> couplings = ising.get_couplings();

    stringstream interaction_stream;
    interaction_stream << "./data_graph/interaction_N" << N << "_T" << fixed << setprecision(2) << T << "_p" << p << (is_dis ? "_dis" : "_ord") << ".dat";
    string interaction_filename = interaction_stream.str();
    ofstream outfile1(interaction_filename);

    stringstream config_stream;
    config_stream << "./data_graph/config_N" << N << "_T" << fixed << setprecision(2) << T << "_p" << p << (is_dis ? "_dis" : "_ord") << ".dat";
    string config_filename = config_stream.str();
    ofstream outfile2(config_filename);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            outfile1 << couplings[i][j] << "\t";
        }
        outfile1 << "\n"; 
    }

    for(int j=0; j<NUM_SWEEPS; j++){
        ising.montecarlo_sweep(); // Sweep di Montecarlo
        vector<int> spins = ising.get_spins();
        for (int k = 0; k < N; k++) {
            outfile2 << spins[k] << "\t"; 
        }
        outfile2 << endl;
    }
    return 129;
    #endif

}