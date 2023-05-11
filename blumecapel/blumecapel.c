/*by Stefano Bae*/
//Simulazione Blume-Capel per reticolo 2D o ER Random Graph
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_matrix.h>

//Modificare i commenti per selezionare il reticolo o il grafo random
#define LAT  // Lattice
//#define ERRG //ER Random Graph

//Modificare THERM per impostare il numero di passi per la termalizzazione iniziale
#define THERM 100000

//Modificare NNRANGE per impostare il numero massimo di vertici
#define NNRANGE 10

//Modificare per attivare l'output in formato .txt
//#define TXTOUT

/*----------Variabili Globali----------*/
double inv3 = 1./3.;
int N, L;           //Numero di spins, dimensione del reticolo
int intmat = 0;     
double J;           //Costante di interazione
double T, beta;     //Temperatura e temperatura inversa
double mu;          //Potenziale chimico
double ***pacc;     //Look up table per le probabilita' di accettazione
double p;           

void exit_failure(char *s){
    fprintf(stderr,"%s",s);
    exit(EXIT_FAILURE);
}

double drand(){
    return (double) rand()/RAND_MAX;
}

void init_rand(){
    //Cambiare con qualcosa di migliore
    srand(time(0));
}

void init_config(int N,int *s, gsl_matrix_int *map,gsl_matrix *Jmat){
    double r;
    #ifdef LAT 
    int n,m;
    for(int i = 0;i<N;i++){
        if(intmat==0){
            //map near neighbours
            n = i%L;
            m = floor(i/L);

            gsl_matrix_int_set(map,i,0,(N+i-1)%N);
            gsl_matrix_int_set(map,i,1,(N+i+1)%N);
            gsl_matrix_int_set(map,i,2,n+L*((L+m-1)%L));
            gsl_matrix_int_set(map,i,3,n+L*((L+m+1)%L));
            gsl_matrix_int_set(map,i,4,-1);

            gsl_matrix_set(Jmat,i,(N+i+L)%N,J/T);
            gsl_matrix_set(Jmat,i,(N+i-1)%N,J/T);
            gsl_matrix_set(Jmat,i,(N+i+1)%N,J/T);
            gsl_matrix_set(Jmat,i,(N+i-L)%N,J/T);
        }
        //assign random value of spin
        r = drand();
        if(r<inv3) s[i] =-1;
        if((r>=inv3)&&(r<2.*inv3)) s[i] = 0;
        if(r>= 2.*inv3) s[i] = 1;
    }
    #endif
    #ifdef ERRG
    int count;
    int countj[N];
    for(int i = 0;i<N;i++) countj[i]=0;

    for(int i = 0;i<N;i++){
        if(intmat == 0){
            count = 0;
            for(int j = i+1;j<N;j++){
                //Given the i-th spin, loop over the other spins
                //If p>r add the link between the i-th and the j-th spin
                r = drand();
                if(p>r){
                    gsl_matrix_int_set(map,i,countj[i]+count,j);
                    gsl_matrix_int_set(map,j,countj[j],i);
                    /*map[i][countj[i]+count] = j;
                    map[j][countj[j]] = i;*/
                    gsl_matrix_set(Jmat,i,j,J/T);
                    gsl_matrix_set(Jmat,j,i,J/T);
                    
                    count++;
                    countj[j]++;
                }
            }
            gsl_matrix_int_set(map,i,countj[i]+count,-1);
        }
        //map[i][countj[i]+count]=-1;
        //assign the i-th spin a random value
        r = drand();
        if(r<inv3) s[i] =-1;
        if((r>=inv3)&&(r<2.*inv3)) s[i] = 0;
        if(r>= 2.*inv3) s[i] = 1;
    }
    #endif
}

int new_spin(int spin){
    return spin + 3*(spin == -2) - 3*(spin == 2);
}

void init_obs(int N,double *energy, int *magn, int *rho, int *s, gsl_matrix_int *map){
    double tmp_e = 0.;
    int tmp_rho,tmp_m;
    int nnsum =0;
    tmp_rho = tmp_m = 0;
    int count=0;

    for(int i = 0;i<N;i++){
        count = 0;
        nnsum = 0;
        while(gsl_matrix_int_get(map,i,count) != -1){
            nnsum += s[gsl_matrix_int_get(map,i,count)];
            count++;
        }
        tmp_e += -0.5*(double)s[i]*J*nnsum + mu*s[i]*s[i];
        tmp_m += s[i];
        tmp_rho += s[i]*s[i];
    }
    (*energy) = tmp_e;
    (*magn) = tmp_m;
    (*rho) = tmp_rho;
}

void init_pacc(){
    int i,j,k;
    for(i = 0;i<3;i++){
        for(j = 0;j<5;j++){
            for(k = 0;k<(2*NNRANGE+1);k++) pacc[i][j][k] = exp(-(double)beta*(-J*(j-2)*(k-NNRANGE) + mu*(i-1)));
        }
    }
}

void one_sweep_heli(int N,int *s,gsl_matrix_int *map, double *dE, int *dM,int *dRho){
    int idx;
    int count;
    int nnsum;
    int ns,os;
    double deltaE;
    int deltaRho;
    double r;

    for(int i = 0;i<N;i++){
        idx = rand()%N;
        count = 0;
        nnsum = 0;
        while(gsl_matrix_int_get(map,idx,count)!=-1){
            nnsum += s[gsl_matrix_int_get(map,idx,count)];
            count++;
        }
        r = drand();
        os = s[idx];
        ns = new_spin(os + (r<0.5) - (r>=0.5));
        deltaRho = ns*ns-os*os;
        deltaE = -(double)J*(ns-os)*nnsum + mu*deltaRho;

        if(deltaE <= 0.){
            s[idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }else if(pacc[deltaRho+1][(ns-os)+2][nnsum+NNRANGE] > drand()){
            s[idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }

    }
}

int main(int argc, char **argv){
    int count;
    int MCS;
    double en,de;
    int magn, rho, dm, drho;
    double avg_magn, avg_en, avg_rho;
    FILE *fpConfig, *fpObs,*fpJ;
    char fnameConfig[50], fnameObs[50],fnameJ[50];
    gsl_matrix_int *config;
    gsl_matrix_int *map;
    gsl_matrix *Jmat;

    //Impostazione parametri di simulazione
    if(argc != 8)
        exit_failure("Inserire gli argomenti nel seguente ordine: L, J, mu, T, p, MCS, count\n");
    
    L = atoi(argv[1]);
    J = atof(argv[2]);
    mu = atof(argv[3]);
    T = atof(argv[4]);
    p = atof(argv[5]);
    beta = 1./T;
    MCS = atoi(argv[6]);
    count = atoi(argv[7]);

    //printf("Blume Capel:\nL = %d\nJ = %f\nmu = %f\nT = %f\n\n",L,J,mu,T);
    N = L*L;
    int s[N];    
    Jmat = gsl_matrix_calloc(N,N);
    map = gsl_matrix_int_calloc(N,N);
    config = gsl_matrix_int_calloc(MCS,N);
    init_rand();

    //Allocazione memoria
    pacc = (double***)malloc(3*sizeof(double));
    for(int i = 0; i<3;i++){
        pacc[i] = (double**)malloc(5*sizeof(double));
        for(int j = 0;j<5;j++){
            pacc[i][j] = (double*)malloc((2*NNRANGE+1)*sizeof(double));
        }
    }
    
    //Inizializzo probabilita' di accettazione
    init_pacc();

    //Leggi interazioni esistenti
    //Altrimenti genera una matrice di interazione
    sprintf(fnameJ,"outputs/bcMap.bin");
    fpJ = fopen(fnameJ,"rb");
    if(fpJ != NULL){
        gsl_matrix_int_fread(fpJ,map);
        intmat = 1;
        //printf("Read bcMap\n");
        fclose(fpJ);
    }
    
    //Inizializzo configurazione random
    init_config(N,s,map,Jmat);
    if(intmat == 0){
        fpJ = fopen(fnameJ,"wb");
        gsl_matrix_int_fwrite(fpJ,map);
        fclose(fpJ);
        sprintf(fnameJ,"outputs/bcJ0.bin");
        fpJ = fopen(fnameJ,"wb");
        gsl_matrix_fwrite(fpJ,Jmat);
        fclose(fpJ);
        //printf("Init conf\n");
    }

    init_obs(N,&en,&magn,&rho,s,map);

    sprintf(fnameConfig,"outputs/bcConfig_n%d.bin",count);
    sprintf(fnameObs,"outputs/bcObs_L%d_%.4lf_%.4lf.txt",L,T,mu);
    fpConfig = fopen(fnameConfig,"wb+");
    fpObs = fopen(fnameObs,"w+");

    //printf("Energy = %.3lf\nMagn = %.3lf\nRho = %.3lf\n\n",(double)en/N,(double)magn/N,(double)rho/N);
    int i;
    avg_en = avg_magn =avg_rho =0.;
    clock_t begin = clock();
    for(i=0;i<THERM;i++){
        de= 0.;
        drho = dm = 0;
        one_sweep_heli(N,s,map,&de,&dm,&drho);
        magn += dm;
        en += de;
        rho += drho;
    }
    //printf("OK\n");
    for(i=0;i<MCS;i++){
        de= 0.;
        drho = dm = 0;
        one_sweep_heli(N,s,map,&de,&dm,&drho);
        magn += dm;
        en += de;
        rho += drho;
        
        avg_en += en;
        avg_magn += magn;
        avg_rho += rho;

        fprintf(fpObs,"%f\t%f\t%f\n",en/N,(double)magn/N,(double)rho/N);
        for(int k = 0;k<N;k++) gsl_matrix_int_set(config,i,k,s[k]);
        
    }
    //printf("OK\n");
    
    clock_t end =clock();
    double runtime = (double)(end-begin)/CLOCKS_PER_SEC;
    //printf("<e> = %f\t<m> = %f\t<rho> = %f\n",avg_en/(N*MCS),(double)avg_magn/(N*MCS),(double)avg_rho/(N*MCS));
    printf("\nRuntime = %.4f [s]\n",runtime);

    gsl_matrix_int_fwrite(fpConfig,config);
    gsl_matrix_int_free(config);
    gsl_matrix_int_free(map);
    gsl_matrix_free(Jmat);

    fclose(fpObs);
    fclose(fpConfig);
    
    return(EXIT_SUCCESS);
}
