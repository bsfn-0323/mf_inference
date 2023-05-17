#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <time.h>

//correlazione 
double corr2p(int MCS,gsl_vector_int *a,gsl_vector_int *b){
    int mean1 = 0,mean2 = 0;
    int corr = 0;
    for(int i = 0;i<MCS;i++){
        corr += gsl_vector_int_get(a,i)*gsl_vector_int_get(b,i);
        mean1 += gsl_vector_int_get(a,i);
        mean2 += gsl_vector_int_get(b,i);
    }
    return (double)corr/MCS - (double)(mean1*mean2)/(MCS*MCS);
}

void get_corr_mat(int MCS,int N, gsl_matrix_int *config,gsl_matrix *corr){
    double c;
    gsl_vector_int *a = gsl_vector_int_alloc(MCS);
    gsl_vector_int *b = gsl_vector_int_alloc(MCS);
    for(int x =0 ;x<N;x++){
        gsl_matrix_int_get_col(a,config,x);
        c = corr2p(MCS,a,a);
        gsl_matrix_set(corr,x,x,c);
        for(int y = x+1;y<N;y++){
            gsl_matrix_int_get_col(a,config,x);
            gsl_matrix_int_get_col(b,config,y);
            c = corr2p(MCS,a,b);
            gsl_matrix_set(corr,x,y,c);
            gsl_matrix_set(corr,y,x,c);
        }
    }
    gsl_vector_int_free(a);
    gsl_vector_int_free(b);
}

double get_gamma(int N, gsl_matrix *J, gsl_matrix *J0,double T){
    double num = 0.;
    double den = 0.;
    double JJ,JJ0;
    for(int i = 0;i<N;i++){
        for(int j = i+1;j<N;j++){
            JJ = -gsl_matrix_get(J,i,j)/T;
            JJ0 = gsl_matrix_get(J0,i,j);
            num += (JJ-JJ0)*(JJ-JJ0);
            den += JJ0*JJ0;
        }
    }
    //printf("num = %f\nden = %f\n\n",num,den);
    return sqrt(num/den);
}

gsl_matrix *inv_mat(int size,gsl_matrix *matrix){
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;
    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p,&s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);


    gsl_permutation_free(p);

    return inv;
}

int main(int argc, char **argv){
    int MCS,N,L,meas,skip;
    double Tmin,Tmax,dT;
    double gamma;
    FILE *fp;
    gsl_matrix_int *conf;
    gsl_matrix *corr;
    gsl_matrix *J,*J0;
    char fnameConfig[50],fnameJ[50],fnameGamma[50];

    if(argc!=7){
        fprintf(stderr,"Inserire i parametri: L, MCS, Tmin, Tmax, dT, meas \n\n");
        exit(EXIT_FAILURE);
    }
    L = atoi(argv[1]);
    MCS = atoi(argv[2]);
    N = L*L;
    Tmin = atof(argv[3]);
    Tmax = atof(argv[4]);
    dT = atof(argv[5]);
    meas = atoi(argv[6]);
    //skip = atoi(argv[7]);

    conf = gsl_matrix_int_alloc(MCS,N);
    corr = gsl_matrix_alloc(N,N);
    J = gsl_matrix_alloc(N,N);
    J0 = gsl_matrix_alloc(N,N);

    sprintf(fnameGamma,"outputs/bcGamma.txt");
    
    clock_t begin = clock();
    //printf("Reading from %d\n",skip);
    for(int count = 0;count<meas;count++){
        //Leggi configurazione
        sprintf(fnameConfig,"outputs/bcConfig_n%d.bin",count);
        fp = fopen(fnameConfig,"rb");
        gsl_matrix_int_fread(fp,conf);
        fclose(fp);
        //printf("Input Files\n\n");

        //Leggi la matrice di interazione originale
        sprintf(fnameJ,"outputs/bcJ0.bin");
        fp = fopen(fnameJ,"rb");
        gsl_matrix_fread(fp,J0);
        fclose(fp);

        //Calcola la matrice di correlazione
        get_corr_mat(MCS,N,conf,corr);
        //printf("Computed Correlation Matrix\n\n");
        
        //La matrice J inferita, nell'approccio MF, è l'inversa della matrice di correlazione
        J = inv_mat(N,corr);
        //printf("Computed Interaction Matrix\n\n");
        for(int i = 0;i<N;i++) gsl_matrix_set(J,i,i,0.);

        //Salva in binario la matrice di interazione
        sprintf(fnameJ,"outputs/bcJ_n%d.bin",count);
        fp=fopen(fnameJ,"wb");
        gsl_matrix_fwrite(fp,J);
        fclose(fp);
        
        //Compute the reconstruction error
        gamma=get_gamma(N,J,J0,Tmin+dT*count);
        //printf("Computed Gamma\n\n");
        fp=fopen(fnameGamma,"a");
        fprintf(fp,"%f\t%f\n",Tmin+dT*count,gamma);
        fclose(fp);
        
        printf("----------------%d---------------\n",count);
    }
    
    clock_t end =clock();
    double runtime = (double)(end-begin)/CLOCKS_PER_SEC;
    printf("\nRuntime = %.4f [s]\n",runtime);

    gsl_matrix_int_free(conf);
    gsl_matrix_free(corr);
    gsl_matrix_free(J);
    gsl_matrix_free(J0);
    return 1;
}
