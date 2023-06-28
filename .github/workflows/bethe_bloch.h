#ifndef STRUTTURE_FUNZIONI
#define STRUTTURE_FUNZIONI


typedef short unsigned suns; //da 0 a 255, sufficiente per numero atomico

typedef enum targetta
 {
  alluminio=0,
  ferro    =1,
  rame     =2
 }targetta;

typedef enum proiet
 {
  alfa   =0,
  muone  =1,
  protone=2
 }proiet;


typedef struct _target 
 {
  suns   num_atom;   
  double massa_atom; // g (mol)^-1
  double meanexen;   // eV
  double densita;    // g/cm^3
  double spessore;   // cm
 } target;  
typedef struct _proiettile
 {
  int    carica;  // unità di e
  double energia; // MeV
  double massa;   // MeV/c^2
 } proiettile;
typedef target*      ptarget;
typedef proiettile*  p_proiettile;


//funzione che dato il target e il suo spessore, alloca in memoria lo spazio per un puntatore a target
//e lo inizializza con i valori dati e quelli tabulati, per riempire la sua struct creata sopra
ptarget give_ptarget(targetta tar, double spessore);

//funzione che dato il proiettile e la sua energia, alloca in memoria lo spazio per un puntatore a proiettile
//e lo inizializza con i valori dati e quelli tabulati, per riempire la sua struct creata sopra
p_proiettile give_p_proiettile(proiet prtt, double energy);

//funzione che, dato proiettile e target, calcola e dunque fornisce in output lo stopping power medio in MeV/cm
double stopping_power_medio(ptarget ptar, p_proiettile p_pr);

//funzione che, dato proiettile, target, e passi di integrazione
//calcola lo stopping power integrale in e fornisce in output l'energia rimanente in MeV
double stopping_power_integrale(ptarget ptar, p_proiettile p_pr, unsigned passi);

//funzione che, dato proiettile, target e passi di integrazione
//calcola lo stopping power integrale in e fornisce in output l'energia rimanente in MeV. 
//Memorizza poi i singoli dE in un file chiamato "dati.txt".
double stopping_power_print(ptarget ptar, p_proiettile p_pr, unsigned passi);


//funzione che, dato proiettile, target, energia finale e dx minimo
//valuta lo spessore necessario del target, procedendo per spessori dx dati dal dx minimo,
//per degradare un fascio di particelle proiettili, fino a portarle all'energia finale indicata.
double perdita_in_range_particella(ptarget ptar, p_proiettile p_pr, double energia_fin, double dx_min);
#endif
