#include <stdio.h>
#include <stdlib.h>
#include <values.h> // Per MAXFLOAT
#include <math.h> //per fer el pow
#define MAXNumlevels 60
void ExitError(const char *lloc, const char *miss, int err) { fprintf (stderr, "\nERROR %s: %s.\nAbortem el procés ...\n\n", lloc, miss); exit(err); }

//-----ESTRUCTURES-----
//1 ESTRUCTURA ARESTA
typedef struct{ 
  unsigned num_node_adjacent;  //num de nodes adjacents
  float temps; //temps o cost
}aresta;

//2 ESTRUCTURA NODE
typedef struct{
  unsigned n_arestes;  //nomnre d'arestes d'un node
  aresta arestes[5];  //
} node;

//3 ESTRUCTURA NODE_ROUTING_STATUS
typedef struct{
  float cost_g; 
  unsigned node_pare; 
} node_routing_status;

//4 ESTRUCTURA CUA PRIO
typedef struct {
  char nlevels;
  unsigned long last_level_size;
  unsigned *Cua_Prio_BiHe_level[MAXNumlevels];
} Cua_Prio_BiHe;
//5 ELEMENT DE LA CUA
typedef struct {
    int level;
    unsigned long index;
} CPBH_Element;

#define NNODES 21

//-----INICIALITZACIÓ DE FUNCIONS-----
int Dijkstra(unsigned, node *, node_routing_status *, unsigned);
int CPBH_EsBuit(Cua_Prio_BiHe *);
unsigned CPBH_desencua(Cua_Prio_BiHe *, node_routing_status *);
int CPBH_encua(unsigned, node_routing_status *, Cua_Prio_BiHe *);
void CPBH_reencua(unsigned, node_routing_status *, Cua_Prio_BiHe *);
void CPBH_heapify_up(CPBH_Element num_node, Cua_Prio_BiHe *BH, node_routing_status *nods);
void CPBH_heapify_down(CPBH_Element num_node, Cua_Prio_BiHe *BH, node_routing_status *nods);
CPBH_Element CPBH_LookUp (unsigned num_node, Cua_Prio_BiHe *BH);
CPBH_Element CPBH_pare(CPBH_Element e, Cua_Prio_BiHe *BH);
CPBH_Element CPBH_right_son(CPBH_Element e, Cua_Prio_BiHe *BH);
CPBH_Element CPBH_left_son(CPBH_Element e, Cua_Prio_BiHe *BH);
float CPBH_getcost(unsigned num, node_routing_status *nods);


//-----MAIN-----
int main() { register unsigned v;
/* Inicialitzacio */
    node_routing_status status_nodes[NNODES];
    node nodes[NNODES] = {
        {3, { {1, 0.528},   {2, 0.495},   {3, 0.471} }},                 // 0
        {2, { {0, 0.528},   {3, 0.508} }},                               // 1
        {4, { {0, 0.495},   {3, 3.437},   {5, 12.033},  {15, 34.852} }}, // 2
        {4, { {0, 0.471},   {1, 0.508},   {2, 3.437},   {4, 23.155} }},  // 3
        {4, { {3, 23.155},  {5, 6.891},   {6, 4.285},   {7, 0.520} }},   // 4
        {2, { {2, 12.033},  {4, 6.8910} }},                              // 5
        {3, { {4, 4.285},   {7, 0.630},   {8, 17.406} }},                // 6
        {2, { {4, 0.520},   {6, 0.630} }},                               // 7
        {5, { {6, 17.406},  {9, 6.657},   {10, 15.216}, {11, 10.625}, {12, 17.320} }}, // 8
        {2, { {8, 6.657},   {12, 16.450} }},                             // 9
        {2, { {8, 15.216},  {14, 12.373} }},                             // 10
        {2, { {8, 10.625},  {12, 3.618} }},                              // 11
        {3, { {8, 17.320},  {9, 16.450},  {11, 3.618} }},                // 12
        {2, { {14, 4.450},  {19, 6.450} }},                              // 13
        {4, { {10, 12.373}, {13, 4.450},  {15, 16.178}, {19, 5.203} }},  // 14
        {5, { {2, 34.852},  {14, 16.178}, {16, 4.818},  {18, 3.877},  {19, 19.131} }}, // 15
        {3, { {15, 4.818},  {17, 3.199},  {18, 2.976} }},                // 16
        {2, { {16, 3.199},  {18, 20.832} }},                             // 17
        {4, { {15, 3.877},  {16, 2.976},  {17, 20.832}, {20, 2.510} }},  // 18
        {4, { {13, 6.450},  {14, 5.203},  {15, 19.131}, {20, 13.313} }}, // 19
        {2, { {18, 2.510},  {19, 13.313} }} };                           // 20
    unsigned node_origen = 0;
    if(Dijkstra(node_origen, nodes, status_nodes, NNODES) == 66) ExitError("a Dijkstra", "no és possible allocatar memòria per la cua", 66);
/* Detecció dels nodes extremals: són els que no són pare de ningú */
    char EsPare[NNODES] = { 0U };
    for (v=1; v < NNODES; v++) {
        if(!(status_nodes[v].cost_g < MAXFLOAT)) ExitError("al graf", "hi ha nodes inaccessibles", 3);
        EsPare[status_nodes[v].node_pare] = 1;
    }

/* Per cada node extremal trobem el camí des-de l'origen */
    unsigned Cami[NNODES];
    for(v=1; v < NNODES; v++) {
        if(EsPare[v]) continue;
        register int l = 0;
        Cami[0] = v;
        do { l++; Cami[l] = status_nodes[Cami[l-1]].node_pare;} while (Cami[l]) ;
        fprintf(stdout, "[0]");
        for(l-- ; l >= 0; l--) fprintf(stdout, " -- %g --> [%d]", status_nodes[Cami[l]].cost_g, Cami[l]) ;
        fprintf(stdout, "\n");
    }
    return 0;
}

/* La fase d'inicialització de Dijkstra inclou establir
 *          nodstatus[u].cost_g := MAXFLOAT;     i     Node_es_a_S[v] = 0;
 * per a tot node v.
 *
 * Per altra banda quan un node entra a la cua (s'encua) ho fa amb nodstatus[v].cost_g < MAXFLOAT
 * i pot fer una de les dues accions següents:
 *     * Moure's dins de la cua (reencua) amb un valor nodstatus[v].cost_g < MAXFLOAT més baix del que tenia fins la iteració actual; o
 *     * Sortir de la cua amb la regla ExtractMin (Desencua). En aquest cas sabem que el seu cost ja és optimal i l'hem d'afegir al conjunt S
 *       En aquest cas tenim nodstatus[u].cost_g < MAXFLOAT i marquem que v ? S establint Node_es_a_S[v] = 1;
 *
 * Llavors, un node arbitrari v pot estar en un dels tres estats següents
 * v ? S (S és el conjunt de nodes dels que ja se sap la distància mínima) <==> Node_es_a_S[v] == 1 (notis que, en aquest cas, tenim nodstatus[u].cost_g < MAXFLOAT)
 * v ? S  i, a més, v no és a la Cua_Prio_BiHe <==> nodstatus[u].cost_g == MAXFLOAT (notis que, en aquest cas, tenim Node_es_a_S[v] == 0)
 * v ? S però v és a la Cua_Prio_BiHe <==> nodstatus[u].cost_g < MAXFLOAT && Node_es_a_S[v] == 0 */
 int Dijkstra(unsigned source, node *nods, node_routing_status *nodstatus, unsigned nnods ){ register unsigned v;
    char Node_es_a_S[nnods];
    for (v=0; v < nnods; v++) { nodstatus[v].cost_g = MAXFLOAT; Node_es_a_S[v] = 0; } // Totes les distàncies a infinit

    Cua_Prio_BiHe LaCua = { 0, 0 };
    nodstatus[source].cost_g = 0.0; nodstatus[source].node_pare = UINT_MAX; // El node source és l'arrel i, per tant, no té pare
    CPBH_encua(source, nodstatus, &LaCua);

// Llaç principal
    while(!CPBH_EsBuit(&LaCua)){ register unsigned e;
        unsigned u = CPBH_desencua(&LaCua, nodstatus); Node_es_a_S[u] = 1;
        for(e=0; e < nods[u].n_arestes; e++){ /* Llaç d'expansió del node u que acabem de posar a S */
            v = nods[u].arestes[e].num_node_adjacent;
            if(Node_es_a_S[v]) continue;
            float cost = nodstatus[u].cost_g + nods[u].arestes[e].temps;
            if(nodstatus[v].cost_g > cost) {
                int IsNodeInLaCua = nodstatus[v].cost_g < MAXFLOAT; /* Relaxation step */
                nodstatus[v].cost_g = cost; nodstatus[v].node_pare = u;
                if(IsNodeInLaCua) CPBH_reencua(v, nodstatus, &LaCua); else { int res_encua = CPBH_encua(v, nodstatus, &LaCua); if(res_encua) return 66; }
            } /* Fi del Relaxation step */
        } /* Fi del llaç d'expansió del node u que acabem de posar a S */
    } /* Fi del Llaç principal */
    return 0;
}

/***********************************
 * Funcions de cua amb Binary Heap *
 ********************************************************************************
 * Recordem la Propietat de forma del BinHeap:
 * el BinHeap és un arbre binari complet; és a dir, tots els nivells de l’arbre,
 * excepte possiblement l’últim (el més profund) estan completament omplerts
 * (amb dos fills per node) i, si l’últim nivell de l’arbre no està complet,
 * els nodes d’aquest nivell s’omplen d’esquerra a dreta.
 *
 * La Propietat de forma del BinHeap traduida a la representació
 *      typedef struct {
 *          char nlevels;
 *          unsigned long last_level_size;
 *          unsigned *Cua_Prio_BiHe_level[MAXNumlevels];
 *      } Cua_Prio_BiHe;
 * s'implmenta d'acord amb la següent
 * CONVENCIÓ: Tots els nivells l (de 0 a nlevels-2) han d'estar plens (amb 2^l elements).
 *            L'últim nivell (l = nlevels-1) ha de ser no buit, però pot no estar complet ==>
 *                    el nombre d'elements del darrer nivell (nlevels-1), denotat per last_level_size,
 *                    ha de verificar 0 < last_level_size <= 2^(nlevels-1)
 *            nlevels = 0 és un heap sense nivells; és a dir, un heap buit
 *            OBSERVEU que, quan el heap és no buit, el nivell 0 sempre està ple, ja que té com a màxim i com a mínim un element.
* MAXNumlevels = 60: Del que s'ha dit abans, tenim que la capacitat màxima total del heap és de 2^(MAXNumlevels)-1 = 2^(60)-1, realment més que suficient per a una cua.
 *
 * A més cal que es compleixi la Propietat Heap:
 * la clau emmagatzemada a cada node és menor o igual que les claus dels fills del node, segons algun ordre total.
 * Iterativament la Propietat heap diu que el cost d'un node donat N és més petit o igual que els valors de tots els nodes del subarbre que té el node N com arrel.
 ******************************************************************************************************************************************************************/

//-----FUNCIONS-----
int CPBH_EsBuit(Cua_Prio_BiHe *BH){ 
  if (BH->nlevels==0){return 1;}//mirem si està buida
  else{
    return 0;
  }
}
//FUNCIÓ DESENCUA

unsigned CPBH_desencua(Cua_Prio_BiHe *BH, node_routing_status *nods){ 
 //per accedir al node ->Cua_prio_element[nivell][node]
  //Per accedir a l'ultim nivell  nlevels (quantitat de nivells que hi ha), i ultim node   last_level_size (cuants nodes hi ha a l'ultim nivell) i últim node l
  //1r pas llegir i retornar el node arrel
  unsigned arrel =BH->Cua_Prio_BiHe_level[0][0];//node arrel
  unsigned ultim_node=BH->Cua_Prio_BiHe_level[BH->nlevels-1][BH->last_level_size-1];//node ultim
  BH->Cua_Prio_BiHe_level[0][0]=ultim_node;//situar l'ultim element a l'arrel
  if(BH->last_level_size>1 ){ //Traiem l'ultim node del nivell, no acaba nivell
    BH->last_level_size=BH->last_level_size-1;//si no acabem amb un nivell, hi ha un node menys a l'ultim mencionat.
  }  
  //cas on triem ultim node del nivell
  else{
    BH->nlevels=BH->nlevels-1; //tenim un nivell menys
    BH->last_level_size=pow(2,BH->nlevels-1); //el nombre d'elements a l'ultim nivell és màxim
  }
  
  if(BH->nlevels!=1){ 
   CPBH_Element el_arrel;
   el_arrel.level=0;
   el_arrel.index=0;
   CPBH_heapify_down(el_arrel,BH,nods);
}
  return arrel; 
  }


//FUNCIO REENCUA
void CPBH_reencua(unsigned num_node, node_routing_status *nods, Cua_Prio_BiHe *BH){
 CPBH_Element element_node =CPBH_LookUp(num_node, BH);//camí més curt al d'abans, es fa un heapify repetidament
    CPBH_heapify_up(element_node, BH, nods);
}

//FUNCIO ENCUA
int CPBH_encua(unsigned num_node, node_routing_status *nods, Cua_Prio_BiHe *BH){
  //Mirar si es buit i si ho és afegirg-ho al [0][0] (AFEGIR EL NODE ARREL)
  if(CPBH_EsBuit(BH)==1){
    BH->nlevels=1;//només hi ha un nivell
    BH->last_level_size=1;//el size de nivell 0 és 1 (arrel)
    BH->Cua_Prio_BiHe_level[(BH->nlevels)-1]=(unsigned*)malloc (pow(2,(BH->nlevels-1))*sizeof(unsigned));
    BH->Cua_Prio_BiHe_level[0][0]=num_node;//afegim aquest num_node al nivell 0 node 0 (l'unic al nivell)
    return 0;//ja hem encuat
  }
  //sino afegir el node al darrer nivell al darrer node consecutiu
  //1r calcular si el nivell esta ple
  else if(BH->last_level_size == pow(2,BH->nlevels-1)){//el nivell que estem mirant és ple
    BH->nlevels=BH->nlevels+1;//Ampliem el nombre de levels que hi ha
    BH->last_level_size = 1;//només hi ha un node al nivell creat (el que acabem d'afegir)
    BH->Cua_Prio_BiHe_level[(BH->nlevels)-1]=(unsigned*)malloc (pow(2,(BH->nlevels-1))*sizeof(unsigned));
    BH->Cua_Prio_BiHe_level[BH->nlevels-1][0]=num_node;//guardem el valor de num_node al seu lloc corresponent
}
  else{
    BH->Cua_Prio_BiHe_level[BH->nlevels-1][BH->last_level_size]=num_node;
    BH->last_level_size=BH->last_level_size+1;
    }
    CPBH_Element el;
    el.level=BH->nlevels-1;
    el.index=BH->last_level_size-1;
    CPBH_heapify_up(el,BH,nods);
    return 0; 
  }


//FUNCIÓ LOOKUP
CPBH_Element CPBH_LookUp (unsigned num_node, Cua_Prio_BiHe *BH){
int trobat=0;
CPBH_Element buscat;
 while(trobat==0){
  for(int nivell=0; nivell<= BH->nlevels; nivell++){//hem de buscar nivell a nivell
    for(int i=0;i<pow(2,nivell);i++){//buscar node a node
      if (BH->Cua_Prio_BiHe_level[nivell][i]==num_node){
        buscat.index=i;//al node trobat li assignem les seves dades
        buscat.level=nivell;
        trobat=1;
        return buscat;
      }
    }
  }
  }
  if (trobat==0){printf("node inexistent");}
   return buscat;
  }
  

//TROBAR EL FILL ESQUERRA D'UN NODE
CPBH_Element CPBH_left_son(CPBH_Element e, Cua_Prio_BiHe *BH){
  int nivell_pare=e.level;//l'element entrat per paràmetre és el pare del fill que volem buscar
  int pos_pare=e.index;
  //comporovar abans que el node e és pare, és a dir que aquest no es trobi a l'ultim nivell
  if(nivell_pare==BH->nlevels-1){
    return e;}
  CPBH_Element fill_esq; //creem un element fill esquerra
  fill_esq.index=pos_pare*2;//si ens dibuixem un arbre amb les posicions es veu molt clar
  fill_esq.level=nivell_pare+1;//un nivell més profund que el pare
  return fill_esq;//retornem l'element fill
}

//TROBAR EL FILL DRET D'UN NODE
CPBH_Element CPBH_right_son(CPBH_Element e, Cua_Prio_BiHe *BH){
  CPBH_Element fill_dret;//creem un element fill dret
  int nivell_pare=e.level;
  int pos_pare=e.index;
  //comporovar abans que el node e és pare, és a dir que aquest no es trobi a l'ultim nivel
  if(nivell_pare==BH->nlevels-1){
    return e;}
    fill_dret.index=(pos_pare*2)+1;//si ens dibuixem un arbre amb les posicions es veu molt clar
    fill_dret.level=nivell_pare+1;//un nivell més profund que el pare
  return fill_dret;
    }
//TROBAR EL PARE D'UN NODE
CPBH_Element CPBH_pare(CPBH_Element e, Cua_Prio_BiHe *BH){
  CPBH_Element node_pare;//creem ekl node pare
  int nivell_fill=e.level;
  int pos_fill=e.index;
  if(nivell_fill==0){//no podem buscar el pare d'el node arrel
    return e;}
  node_pare.index=(pos_fill/2);//si ens dibuixem un arbre amb les posicions es veu molt clar
  node_pare.level=nivell_fill-1;//pare un nivell més alt que el fill
  return node_pare;
  }

float CPBH_getcost(unsigned num_node, node_routing_status *nods){
  return nods[num_node].cost_g;
}

void CPBH_heapify_up (CPBH_Element q, Cua_Prio_BiHe *BH, node_routing_status *nods){
  int acabat=0;
    while(acabat==0){
     int nivell_fill=q.level;
     int pos_fill=q.index;
     CPBH_Element info_pare=CPBH_pare(q,BH);
     unsigned node_pare=BH->Cua_Prio_BiHe_level[info_pare.level][info_pare.index];
     unsigned node_fill=BH->Cua_Prio_BiHe_level[nivell_fill][pos_fill];
      
     double cost_pare=nods[node_pare].cost_g;
     double cost_fill=nods[node_fill].cost_g;
     if (cost_pare>cost_fill){//si el pare té menys cost que el fill, hem de intercanviar-los
        BH->Cua_Prio_BiHe_level[info_pare.level][info_pare.index]=node_fill;//node a la pos del pare és el fill
        BH->Cua_Prio_BiHe_level[nivell_fill][pos_fill]=node_pare;//node  a la posiciço del fill és el pare
        q=info_pare;
       }
     else{
        acabat=69;
      }
    }
  }


void CPBH_heapify_down(CPBH_Element q, Cua_Prio_BiHe *BH, node_routing_status *nods){
  int acabat=0;
  if(CPBH_EsBuit(BH)==1){acabat=1;} //cas arbre buit
  while(acabat==0){
    CPBH_Element  info_fill_esq=CPBH_left_son(q,BH);//busquem el fill dret
    CPBH_Element info_fill_dret=CPBH_right_son(q,BH);//busquem el fill esquerra
    if(info_fill_dret.index==info_fill_esq.index){
      break;
    }
    
    unsigned fill_dret=BH->Cua_Prio_BiHe_level[info_fill_dret.level][info_fill_dret.index];//nodes fills i pare
    unsigned fill_esq=BH->Cua_Prio_BiHe_level[info_fill_dret.level][info_fill_dret.index];
    unsigned pare=BH->Cua_Prio_BiHe_level[q.level][q.index];
    float cost_fill_esq=CPBH_getcost(fill_esq,nods);//calculs del cost
    float cost_fill_dret=CPBH_getcost(fill_dret,nods);
    float cost_pare=CPBH_getcost(pare,nods);
    if(cost_pare>cost_fill_esq && cost_pare>cost_fill_esq){//pare més gran que els dos fills
      if(cost_fill_dret<cost_fill_esq){ //pare més gran que els dos fills pero el fill dret més petit
        BH->Cua_Prio_BiHe_level[q.level][q.index]=fill_dret;
        BH->Cua_Prio_BiHe_level[info_fill_esq.level][info_fill_esq.index]=pare;
        q=info_fill_dret;
        }
      else{//pare més gran que els dos fills pero el fill dret més petit
        BH->Cua_Prio_BiHe_level[q.level][q.index]=fill_esq;
        BH->Cua_Prio_BiHe_level[info_fill_dret.level][info_fill_dret.index]=pare;
        q=info_fill_esq;
        }
      }
     else if(cost_pare>cost_fill_dret){//dret més petit que pare pero esquerra més gran
        BH->Cua_Prio_BiHe_level[q.level][q.index]=fill_dret;
        BH->Cua_Prio_BiHe_level[info_fill_dret.level][info_fill_dret.index]=pare;
        q=info_fill_dret;
      }
        
    else if(cost_pare>cost_fill_esq){//esquerra més petit que pare pero dret més gran
      BH->Cua_Prio_BiHe_level[q.level][q.index]= fill_esq;
      BH->Cua_Prio_BiHe_level[info_fill_esq.level][info_fill_esq.index]=pare;
      q=info_fill_esq;
    }
    
    else{acabat=69;}
}
  }

