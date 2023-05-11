#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*----------Variabili Globali----------*/
double inv3 = 1./3.;
int N, L;
double J;
double T, beta;
double mu;
double ***pacc;


/*----------Funzioni-------------*/

double drand(){
    return (double) rand()/RAND_MAX;
}

void init_rand(){
    //Cambiare con qualcosa di migliore
    srand(time(0));
}

int new_spin(int spin){
    return spin + 3*(spin == -2) - 3*(spin == 2);
}

void init_config(int fill, int *s){
    int i;
    double r;
    if(abs(fill) > 1){
        for(i = 0;i<N;i++){
            r = drand();
            if(r<inv3) s[i] =-1;
            if((r>=inv3)&&(r<2.*inv3)) s[i] = 0;
            if(r>= 2.*inv3) s[i] = 1;
        }
    }else{
        for(i = 0;i<N;i++) s[i] = fill;
    }    

}

void print_config(int wfile, int *s){
    int i;
    for(i = 0;i<N;i++){
        if(i%L == 0) printf("\n");
        printf("%d ",s[i]);
    }
    printf("\n\n");
}

void init_obs(double *energy, int *magn, int *rho, int *s){
    double tmp_e = 0.;
    int tmp_rho,tmp_m;
    int nnsum,n,m;
    tmp_rho = tmp_m = 0;
    int i;
    for(i = 0;i<N;i++){
        m = floor(i/L);
        n = i%L;
        nnsum = s[(N+i-1)%N] + s[(N+i+1)%N] + s[n+L*((L+m-1)%L)] + s[n+L*((L+m+1)%L)];
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
            for(k = 0;k<9;k++) pacc[i][j][k] = exp(-(double)beta*(-J*(j-2)*(k-4) + mu*(i-1)));
        }
    }
}

void one_sweep_heli(int *s, double *dE, int *dM,int *dRho){
    int idx;
    int n,m,i;
    int nnsum;
    int ns,os;
    double deltaE;
    int deltaRho;
    double r;

    for(i = 0;i<N;i++){
        idx = rand()%N;
        n = idx%L;
        m = floor(idx/L);

        nnsum = s[(N+idx-1)%N] + s[(N+idx+1)%N] + s[n+L*((L+m-1)%L)] + s[n+L*((L+m+1)%L)];

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
        }else if(pacc[deltaRho+1][(ns-os)+2][nnsum+4] > drand()){
            s[idx] = ns;
            *dM += ns-os;
            *dE += deltaE;
            *dRho += deltaRho;
        }

    }
}

int main(int argc, char *argv[]){
    int fill;
    int MCS;
    double en,de;
    int magn, rho, dm, drho;
    double avg_magn, avg_en, avg_rho;
    double real_pacc[3][5][9];
    FILE *fpConfig, *fpObs;
    char fnameConfig[50], fnameObs[50];
    //Parametri iniziali
    //{L,J,mu} {T} {fill} {MCS}
    if(argc != 7){
        fprintf(stderr, "Error\n");
        exit(EXIT_FAILURE);
    }
    L = atoi(argv[1]);
    J = atof(argv[2]);
    mu = atof(argv[3]);
    T = atof(argv[4]);
    beta = 1./T;
    fill = atoi(argv[5]);
    MCS = atoi(argv[6]);

    printf("Blume Capel 2D:\nL = %d\nJ = %f\nmu = %f\nT = %f\n\n",L,J,mu,T);

    N = L*L;
    int s[N];    
    /*******************************/
    init_rand();
    pacc = (double***)malloc(3*sizeof(double));
    for(int i = 0; i<3;i++){
        pacc[i] = (double**)malloc(5*sizeof(double));
        for(int j = 0;j<5;j++){
            pacc[i][j] = (double*)malloc(9*sizeof(double));
            
        }
    }
    
    init_pacc();
    init_config(fill,s);
    //print_config(0,s);
    init_obs(&en,&magn,&rho,s);

    sprintf(fnameConfig,"bcConfig_L%d_%.4lf_%.4lf.txt",L,T,mu);
    sprintf(fnameObs,"bcObs_L%d_%.4lf_%.4lf.txt",L,T,mu);
    fpConfig = fopen(fnameConfig,"w");
    fpObs = fopen(fnameObs,"w");

    printf("Energy = %.3lf\nMagn = %.3lf\nRho = %.3lf\n\n",(double)en/N,(double)magn/N,(double)rho/N);

    int i;
    avg_en = avg_magn =avg_rho =0.;
    clock_t begin = clock();
    for(i=0;i<MCS;i++){
        de= 0.;
        drho = dm = 0;
        one_sweep_heli(s,&de,&dm,&drho);
        magn += dm;
        en += de;
        rho += drho;
        if(i>(int)0.75*MCS){
            avg_en += en;
            avg_magn += magn;
            avg_rho += rho;
            fprintf(fpObs,"%f\t%f\t%f\n",en/N,(double)magn/N,(double)rho/N);
            for(int k = 0;k<N;k++) fprintf(fpConfig,"%d\t",s[k]);
            fprintf(fpConfig,"\n");
        }
    }
    
    clock_t end =clock();
    double runtime = (double)(end-begin)/CLOCKS_PER_SEC;
    printf("<e> = %f\t<m> = %f\t<rho> = %f\n",avg_en/(N*MCS),(double)avg_magn/(N*MCS),(double)avg_rho/(N*MCS));
    printf("\nRuntime = %.4f [s]\n",runtime);

    fclose(fpObs);
    fclose(fpConfig);
    
    return(EXIT_SUCCESS);
}