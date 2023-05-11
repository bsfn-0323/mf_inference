#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Parametri utilizzati nella simulazione
//L larghezza del reticolo, SIZE numero di siti (pari a LxLxL), SQUARE taglia del quadrato, J parametro
// NUM_SWEEPS è il numero di sweep montecarlo (per Metropolis 1 sweep corrisponde a SIZE passi)
#define L 3
#define TERMALIZE_SWEEPS 100000
#define NUM_SWEEPS 100000

//Define per la interaction matrix
#undef REGULAR_LATTICE
#define ER_GNM_GRAPH

//interaction: matrice delle interazioni J_ij
//neighbours: matrice dei vicini 

int interaction_matrix[L*L*L][L*L*L];
int neighbours_matrix[L*L*L][L*L*L];

const int SIZE =  L*L*L;
const int SQUARE = L*L;

//lookup table necessaria per i prodotti scalari
int const look_up[7] = {0, -1 , 0, 1 , 0, -1, 0};
//X e Y: coordinate degli spin con orientazione n
int const X[] = {1,0,-1,0};
int const Y[] = {0,1,0,-1}; 

//DEFINE per le varie opzioni di simulazione  
#define START_AT_INFINITE_TEMPERATURE
#undef START_AT_ZERO_TEMPERATURE


//struct che contiene le grandezze del sistema cui siamo interessati
//Per ora considero soltanto la magnetizzazione, quindi una strct non è necessaria, però può servire se siamo interessati ad altre grandezze, come ad esempio l'energia
//n[i] è il numero di spin con colore i 
struct quantities {
  int n[4];
};

void metropolis_sweep(int lattice[], double coefficients[], struct quantities* grandezze);
void zero_temperature_setup(int lattice[], int start_spin);
void infinite_temperature_setup(int lattice[], struct quantities* grandezze);

void lattice_zero(void);
void lattice_3D(void);
void lattice_GNM_ER(void);

//###################################################################

int main (int argc, char *argv[]){
  int i,j,  lattice[L*L*L];
  double temperature, J, coefficients[13];
  struct quantities grandezze;
  char file_name[100];
  FILE* fp;

//Setup del generatore di numeri casuali
//Se RG utilizziamo un seed fissato per avere consistency sulle diverse run
#ifdef REGULAR_LATTICE
  unsigned int myseed;
  FILE *devran = fopen("/dev/random","r");
  //FILE *devran = fopen("/dev/urandom","r"); // lower quality
  fread(&myseed, 4, 1, devran);
  fclose(devran);
  srand48(myseed);
#else
  srand48(29);
#endif

  if (argc != 3) {
    fprintf(stderr, "usage: %s <T> <J>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  temperature = atof(argv[1]);
  J = atof(argv[2]);
  
  for(i = 0; i <= 12; i++){
    //Look up table per Metropolis
    //Lo swing maggiore che possiamo avere è di 12J, ad esempio se il centrale è parallelo ai sei adiacenti e viene flippato di 180 gradi
    //In teoria il termine 0 della somma non si dovrebbe usare, è più comodo lasciarlo per individuare gli indici in seguito.
    //ho messo anche i coefficienti dispari, oltre a quelli pari, perché così l'algoritmo è facilmente modificabile per il modello di Potts non vettoriale (nella lookup table definita in precedenza basta effettuare la sostituzione -1 -> 0)
    coefficients[i] = exp(-i*J/temperature);
  }

#ifdef START_AT_ZERO_TEMPERATURE
  //Se partiamo dalla temperatura 0, impostiamo tutti gli spin inizialmente con il colore 0
  zero_temperature_setup(lattice, 0);
  grandezze.n[0] = SIZE;
  grandezze.n[1] = 0;
  grandezze.n[2] = 0;
  grandezze.n[3] = 0;
#endif
  
#ifdef START_AT_INFINITE_TEMPERATURE
  //Se partiamo da temperatura infinita, orientiamo casualmente gli spin
  infinite_temperature_setup(lattice, &grandezze);
#endif

#ifdef REGULAR_LATTICE
  lattice_3D();
#elif defined ER_GNM_GRAPH
  lattice_GNM_ER();
#endif
//Salviamo su file le matrici di vicini e interazione
  sprintf(file_name, "./data/interaction_L%d_T%.2f_J%.2f.dat", L, temperature, J);
  if((fp = fopen(file_name, "w+")) == NULL){
    printf("Errore nell'apertura del file!\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "# Data used: L = %d, T = %.2f, J = %.2f\n\n", L, temperature, J);
  for(i=0; i < SIZE; i++){
    for(j = 0; j < SIZE; j++){
      fprintf(fp, "%d ", interaction_matrix[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  sprintf(file_name, "./data/neighbours_L%d_T%.2f_J%.2f.dat", L, temperature, J);
  if((fp = fopen(file_name, "w+")) == NULL){
    printf("Errore nell'apertura del file!\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "# Data used: L = %d, T = %.2f, J = %.2f\n\n", L, temperature, J);
  for(i=0; i < SIZE; i++){
    for(j = 0; j < SIZE; j++){
      if(neighbours_matrix[i][j] == -1){
        break;
      }
      fprintf(fp, "%d ", neighbours_matrix[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
//Apriamo il file delle configurazioni
  sprintf(file_name, "./data/config_L%d_T%.2f_J%.2f.dat", L, temperature, J);
  //if((fp = fopen(file_name, "wb")) == NULL){   //Alternativa per scrivere in binario
  if((fp = fopen(file_name, "w+")) == NULL){
    printf("Errore nell'apertura del file!\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp, "# Data used: L = %d, T = %.2f, J = %.2f\n#Iter site value\n\n", L, temperature, J);

//Termalizzazione
  for(i = 1; i <= TERMALIZE_SWEEPS; i++){
    metropolis_sweep(lattice, coefficients, &grandezze); 
  }

//Presa dati
  for(i = 0; i < NUM_SWEEPS; i++){
    metropolis_sweep(lattice, coefficients, &grandezze); 
      for(j = 0; j < SIZE; j++){
	      fprintf(fp,"%d %d %d %d\n", i, j, X[lattice[j]], Y[lattice[j]]);
    }
  }
  fclose(fp);
  return 0;
}

//######################################################################



void zero_temperature_setup(int lattice[], int start_spin){
  int site;
  for(site = 0; site < SIZE; site++){
    lattice[site] = start_spin;
  }
}

void infinite_temperature_setup(int lattice[], struct quantities* grandezze){
  int site;
  double rand_num;

  grandezze->n[0]= 0;
  grandezze->n[1] = 0;
  grandezze->n[2] = 0;
  grandezze->n[3] = 0;
  for(site = 0; site < SIZE; site++){
    rand_num = drand48();
    if (rand_num < 0.25){
      lattice[site] = 0;
      grandezze->n[0]++;
    } else if (rand_num < 0.5) {
      lattice[site] = 1;
      grandezze->n[1]++;
    } else if (rand_num < 0.75){
      lattice[site] = 2;
      grandezze->n[2]++;
    } else {
      lattice[site] = 3;
      grandezze->n[3]++;
    }
  }
}

void metropolis_sweep(int lattice[], double coefficients[], struct quantities* grandezze){
  int step, selected_site, delta, old, new, next, i;
  for(step = 0; step < SIZE; step++){
    selected_site = (int)floor(drand48()*SIZE);
    old = lattice[selected_site]; 
    new = (lattice[selected_site] + (int)floor(drand48()*3) + 1) % 4;
    delta = 0;
    next = neighbours_matrix[selected_site][0];
    i = 1;
    while (next >= 0)
    { 
      delta += look_up[lattice[next] - old + 3] - look_up[lattice[next] - new + 3];
      next = neighbours_matrix[selected_site][i++];
    }
    
    if (delta <= 0){
      lattice[selected_site] = new;
      grandezze->n[new]++;
      grandezze->n[old]--;
    }else{
      if(drand48() <= coefficients[delta]){
        lattice[selected_site] = new;
        grandezze->n[new]++;
        grandezze->n[old]--;
      }
    }
  }
}

void lattice_zero(void){
  int i, j;
  for(i = 0; i < SIZE; i++){
    for(j = 0; j < SIZE; j++){
      interaction_matrix[i][j] = 0;
      neighbours_matrix[i][j] = -1;
    }
  }
}

void lattice_3D(void){
  int i,j,k;
  int current, count;
  
  lattice_zero();

  for(i = 0; i < L; i++){
    for(j = 0; j < L; j++){
      for(k = 0; k < L; k++){
        current = SQUARE*i + L*j + k; 

        interaction_matrix[current][SQUARE*i + L*j + (k+1)%L] = 1;
        neighbours_matrix[current][0] = SQUARE*i + L*j + (k+1)%L;
        count = 1;

        if (!interaction_matrix[current][SQUARE*i + L*j + (k-1+L)%L]){
          interaction_matrix[current][SQUARE*i + L*j + (k-1+L)%L] = 1;
          neighbours_matrix[current][count] = SQUARE*i + L*j + (k-1+L)%L;
          count++;
        }

        if (!interaction_matrix[current][SQUARE*i + L*((j+1)%L) + k]){
          interaction_matrix[current][SQUARE*i + L*((j+1)%L) + k] = 1;
          neighbours_matrix[current][count] = SQUARE*i + L*((j+1)%L) + k;
          count++;
        }

        if (!interaction_matrix[current][SQUARE*i + L*((j-1+L)%L) + k]){
          interaction_matrix[current][SQUARE*i + L*((j-1+L)%L) + k] = 1;
          neighbours_matrix[current][count] = SQUARE*i + L*((j-1+L)%L) + k;
          count++;
        }

        if (!interaction_matrix[current][SQUARE*((i+1)%L) + L*j + k]){
          interaction_matrix[current][SQUARE*((i+1)%L) + L*j + k] = 1;
          neighbours_matrix[current][count] = SQUARE*((i+1)%L) + L*j + k;
          count++;
        }

        if (!interaction_matrix[current][SQUARE*((i-1+L)%L) + L*j + k]){
          interaction_matrix[current][SQUARE*((i-1+L)%L) + L*j + k] = 1;
          neighbours_matrix[current][count] = SQUARE*((i-1+L)%L) + L*j + k;
          count++;
        } 
      }
    }
  }
}

void lattice_GNM_ER(void){
  int couples = 3*SIZE;
  int i, selected, selected2, j;

  lattice_zero();
    
  for(i = 0; i < couples; i++){
    do
    {
      selected = (int)floor(drand48()*SIZE);
      do{
        selected2 = (int)floor(drand48()*SIZE);
      } while(selected == selected2);
    } while (interaction_matrix[selected][selected2] != 0);
    interaction_matrix[selected][selected2] = 1;
    interaction_matrix[selected2][selected] = 1;
    j = 0;
    while(neighbours_matrix[selected][j++] >= 0);
    neighbours_matrix[selected][j-1] = selected2;
    j = 0;
    while(neighbours_matrix[selected2][j++] >= 0);
    neighbours_matrix[selected2][j-1] = selected;
  }
}
