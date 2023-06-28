#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"./bethe_bloch.h"
#include"./my_defines.h"


typedef FILE*        pfile;
typedef ptarget      (*tipo_give_ptarget)(double spe);     
typedef p_proiettile (*tipo_give_p_proiettile)(double en);



//IMPLEMENTAZIONI STANDARD
double stopping_power_medio(ptarget ptar, p_proiettile p_pr)
 {
  unsigned z;
  double   gamma, beta, wmax, i, wmi, sp;

  z=p_pr->carica;
  z*=z;
  gamma=(p_pr->energia)/((p_pr)->massa)+1;
  gamma*=gamma;
  beta=1-(1/gamma);              
  wmax=2*MASSA_BETA*beta*gamma;
  i=(ptar->meanexen)*1e-6;
  wmi=wmax/i;
  wmi*=wmi;
  //attento: da qua in poi z, wmax, gamma, beta, i sono già al quadrato
  sp=(ptar->densita)*K_STOPPING_P*z*(ptar->num_atom)/(ptar->massa_atom)*((0.5*log(wmi)/beta)-1);  

  return sp;
 }

double stopping_power_integrale(ptarget ptar, p_proiettile p_pr, unsigned passi)
 {
  unsigned z, j;
  double   gamma, beta, wmax, i, wmi, costante_int, dx, energy, massa_pr,sp;

  massa_pr=p_pr->massa;
  i=(ptar->meanexen)*1e-6;
  z=p_pr->carica;
  z*=z;
  costante_int=(ptar->densita)*K_STOPPING_P*z*(ptar->num_atom)/(ptar->massa_atom);

  for(j=0, dx=(ptar->spessore)/passi, energy=(p_pr->energia); j<passi; j++)
   {
    gamma=(energy/massa_pr)+1;
    gamma*=gamma;
    beta=1-(1/gamma);   
    wmax=2*MASSA_BETA*beta*gamma;
    wmi=wmax/i;
    wmi*=wmi;
    //attento: da qua in poi z, wmax, gamma, beta, i sono già al quadrato
    if( ( sp=fabs(costante_int*((0.5*(log(wmi))/beta)-1)) ) < SP_LIMIT ) {energy=0; break;}   
    energy-=(sp*dx);
   }
  return energy;
 }

double stopping_power_print(ptarget ptar, p_proiettile p_pr, unsigned passi)
 {
  unsigned z, j;
  double   gamma, beta, wmax, i, wmi, costante_int, dx, energy, massa_pr, sp;
  pfile    pf;

  if( (pf=fopen("dati.txt", "w"))==NULL )  {printf("errore in apertura file"); return 0;}
  massa_pr=p_pr->massa;
  i=(ptar->meanexen)*1e-6;
  z=p_pr->carica;
  z*=z;
  costante_int=(ptar->densita)*K_STOPPING_P*z*(ptar->num_atom)/(ptar->massa_atom);
  
  for(j=0, dx=(ptar->spessore)/passi, energy=(p_pr->energia); j<passi; j++)
   {
    gamma=(energy/massa_pr)+1;
    gamma*=gamma;
    beta=1-(1/gamma); 
    wmax=2*MASSA_BETA*beta*gamma;
    wmi=wmax/i;
    wmi*=wmi;
    //attento: da qua in poi z, wmax, gamma, beta, i sono già al quadrato
    sp=fabs( costante_int*((0.5*(log(wmi))/beta)-1) );
    energy-=(sp*dx);
    if ( (sp)<SP_LIMIT ) 
     {
      fprintf(pf,"%lf\n", energy);
      energy=0;
      break;
     } 
    fprintf(pf,"%lf\n", energy );
   }
  for(;j<passi;j++) {fprintf(pf,"0.000000\n");}
  fclose(pf);
  return energy;
  
 }


double perdita_in_range_particella(ptarget ptar, p_proiettile p_pr, double energia_fin, double dx_min)
 {
  unsigned z, j;
  double   gamma, beta, wmax, i, wmi, costante_int, energy, massa_pr,sp;

  massa_pr=p_pr->massa;
  i=(ptar->meanexen)*1e-6;
  z=p_pr->carica;
  z*=z;
  costante_int=(ptar->densita)*K_STOPPING_P*z*(ptar->num_atom)/(ptar->massa_atom);

  for(j=0, energy=(p_pr->energia); ; j++)
   {
    gamma=(energy/massa_pr)+1;
    gamma*=gamma;
    beta=1-(1/gamma);   
    wmax=2*MASSA_BETA*beta*gamma;
    wmi=wmax/i;
    wmi*=wmi;
    //attento: da qua in poi z, wmax, gamma, beta, i sono già al quadrato
    sp=fabs(costante_int*((0.5*(log(wmi))/beta)-1));
    if( (energy-=(sp*dx_min))<energia_fin ) break;
   }
  return j*dx_min;
  }


//IMPLEMENTAZIONI CON PUNTATORI
//target
ptarget give_ptarget_al(double spess)
 {
  ptarget      palluminio;
  
  if( (palluminio=(ptarget)malloc(sizeof(target)))==NULL) return NULL;
  palluminio->num_atom=Z_AL;   
  palluminio->massa_atom=A_AL; 
  palluminio->meanexen=I_AL;   
  palluminio->densita=P_AL;    
  palluminio->spessore=spess; 

  return palluminio;
 }

ptarget give_ptarget_fe(double spess)
 {
  ptarget      pferro;
  
  if( (pferro= (ptarget) malloc(sizeof(target)))==NULL) return NULL;
  pferro->num_atom=Z_FE;   
  pferro->massa_atom=A_FE; 
  pferro->meanexen=I_FE;   
  pferro->densita=P_FE;    
  pferro->spessore=spess; 

  return pferro;
 }

ptarget give_ptarget_cu(double spess)
 {
  ptarget      prame;
  
  if( (prame= (ptarget) malloc(sizeof(target)))==NULL) return NULL;
  prame->num_atom=Z_CU;   
  prame->massa_atom=A_CU; 
  prame->meanexen=I_CU;   
  prame->densita=P_CU;    
  prame->spessore=spess; 

  return prame;
 }

//proiettile
p_proiettile give_p_proiettile_alfa(double energy)
 {
  p_proiettile palfa;
  
  if( (palfa=(p_proiettile)malloc(sizeof(proiettile)))==NULL) return NULL;
  palfa->carica=CARICA_ALFA;
  palfa->massa = MASSA_ALFA;
  palfa->energia=energy;

  return palfa;
 }

p_proiettile give_p_proiettile_muone(double energy)
 {
  p_proiettile pmuon;
  
  if( (pmuon=(p_proiettile)malloc(sizeof(proiettile)))==NULL) return NULL;
  pmuon->carica=CARICA_MUONE;
  pmuon->massa = MASSA_MUONE;
  pmuon->energia=energy;

  return pmuon;
 }

p_proiettile give_p_proiettile_protone(double energy)
 {
  p_proiettile p_protone;
  
  if( (p_protone=(p_proiettile)malloc(sizeof(proiettile)))==NULL) return NULL;
  p_protone->carica=CARICA_PROT;
  p_protone->massa = MASSA_PROT;
  p_protone->energia=energy;

  return p_protone;
 }

tipo_give_p_proiettile give_p_proiettile_arr[]={give_p_proiettile_alfa,give_p_proiettile_muone,give_p_proiettile_protone};
tipo_give_ptarget      give_ptarget_arr[]={give_ptarget_al, give_ptarget_fe, give_ptarget_cu};


//FUNZIONI IN BETHE_BLOCH.H
ptarget give_ptarget(targetta tar, double spe)
 {return give_ptarget_arr[tar](spe);}
p_proiettile give_p_proiettile(proiet prtt, double energy)
 {return give_p_proiettile_arr[prtt](energy);}
