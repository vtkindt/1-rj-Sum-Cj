//#include "Stat.h"
#include <assert.h>
#include <ctype.h>
#include <sys/types.h>
#include <Windows.h>
#include "psapi.h"
#include <stdio.h>
#include <time.h>
#include <process.h>
#include "heapsort.h"
#include "gerarb.h"
//#include "database.h"
//! Remember to update macros in database.h
#define DBMaximumDimension 5500000	//! initial 200000. 2M is reached on 150 jobs. Can set 6M
#define max(x,y) ((x>y)?(x):(y))
#define min(x,y) ((x<y)?(x):(y))
#define R(t,i)	max(t,ed[i])	//! Real Release time
#define F(t,i)  (R(t,i)+pt[i])	//! Real Finishing time
#define NB	250		//! nb jobs, initial 500. Better put 200
#define NB_MAX	NB //! max nb jobs, initial 500. Better put 200
#define	NOI	31		//! nb bit per word
#define NBH	6		//! not used
#include <memory.h>
#define NBBIT	(NB_MAX/NOI+((NB_MAX%NOI>0)?1:0))	//! nb word for NB_MAX bit manip. Only used by stat
typedef struct hebi{
	int niv, brn, binf, temp, *tache, stat[NBBIT];	//! tache contains job id. Later sorted to ert order. First niv positions are dp prefix
	struct hebi *next;								//! brn is the Csum of the fixed part.
} PTR;
PTR	*premier[NB_MAX], *buf, *last = NULL;			//! premier[n] is possibly the start of list of nodes of size n. But only premier[0] is used in depth first, as the beginning of the whole list
int	bsup, arb[NB], qu[NB], arb1[NB], deb1[NB], prov1[NB], rsc; //! rsc should be the current best sol
char	*malloc();
int	nbatt, nbatt1, nbsepa, nbpas, somrel, C1[NB_MAX], nbPolluted;
long	clock(), topt;
short	etat[NB_MAX];
unsigned nb_task;
int 	ed[NB], pt[NB];	//! release date et processing time
unsigned int LBroot;

//==================================
//! Helper functions. Added for Memo
//==================================
#define INT_MAX 2147483647
#define LOG(...)			//printf(__VA_ARGS__)
#define LOG_FOR(i,n,...)	//{for(int i=0; i<n; ++i)printf(__VA_ARGS__);}
#define LOG_IF(cond,...)	//{if(cond)printf(__VA_ARGS__);}
#define PAUSE_ON(cond)		//{if(cond){printf("%d : "#cond"! Pause...\n", __LINE__); getchar();}}
#define FOR_E(i,n)			//for(int i=0; i<n; ++i)
#define LOG_NODE(node)		//{if(node==NULL)LOG("NULL\n"); else \
							//{LOG(#node" : t0=%d;\tniv=%d;\tbrn=%d;\tlb=%d;\trsc=%d\nJobs : \n", node->temp, node->niv, node->brn, node->binf, rsc);\
							//LOG_FOR(i29, node->niv, "%d\t", node->tache[i29]); LOG("#\t");\
							//LOG_FOR(i29, nb_task-node->niv, "%d\t", node->tache[node->niv+i29]);}}
void logNode(PTR * node){
	int i = 0;
	printf("t0=%d;\tniv=%d;\tbrn=%d;\tlb=%d;\trsc=%d\nJobs : \n", node->temp, node->niv, node->brn, node->binf, rsc); 
	for (;i<node->niv;++i)
		printf("(%d,%d,%d)\t", node->tache[i], ed[node->tache[i]], pt[node->tache[i]]);
}
//! best[l] should always be updated as the currently best solution of the current node on level l. 
//  For a non-leaf node, this value is updated by children nodes on level l+1, for leaf nodes it is updated by the node itself.
//! it should ONLY and ALWAYS be updated when ub is improved by a solution. Dom conditions must be disabled in order to not loss the opt node
int best[NB + 1];//!血与泪 = { INT_MAX };
//! used like best but for storing the best lb of a node, computed as the least of all lb of its children. 
//  Note that the lb of the node itself should not be counted if it has children, since it is the lowest.
//! All ending nodes must be considered, any miss would lead to a too high lb so that the opt sol may be missed. Dom conditions must be disabled.
int bestlb[NB + 1];
char memlb[NB + 1] = {0};  //! flag indicating for the current node, we should memo the opt sol or lb. Lb is considered when: a polluting global dominance condition is applied.
BOOL Stra6LBMemoOn = TRUE; //! whether enable lb memo for strategy 6
//! when a node is solved, update it's father's best value. Save the current node if not saved before
//! PREREQUEST: best[niv] must store the optimal value of the current node
//! First piege found: the UB rsc when first computed before branching, may already have a good value. On a node, its children may either be cut or not cut but solved directly with a value worse than rsc. 
//		In this case the optimal value of THAT node cannot be obtained since its optimal node may be in the cut subtree. Basicly the whole node is implicitely cut.
//! When this is called, must be sure that what we got is the OPTIMAL sol
void UpdateBestOnSolve(PTR* node, BOOL saveToDb, BOOL isLb){
	if (memlb[node->niv] == 1)isLb = TRUE;
	if (node->niv < 1 || (isLb && !Stra6LBMemoOn)) return;
	int sol =  best[node->niv] ;//node could be a sentinel
	if (isLb){
		assert((sol == INT_MAX || memlb[node->niv]==1) && bestlb[node->niv]<INT_MAX);
		sol = bestlb[node->niv];
	}
	else assert(sol != INT_MAX);
	int Ci = 0;//! Cost value of the last prefixed job
	for (int i = 0; i <= node->niv-1; ++i){
		Ci = max(Ci, ed[node->tache[i]]) + pt[node->tache[i]];
	}
	//PAUSE_ON(Csumi[node->niv-1] != node->brn);
	//update father's best
	int i = node->niv-1;			
	int t0 = max(Ci, ed[node->tache[i + 1]]);	//! t0 of the remaining part
	if (saveToDb ){
		i = node->niv;
		//check sol (not considering lb memo)
		if (FALSE && !isLb)//node->niv == 33)// && sol==28823)
		{
			FILE*f = fopen("d.dat", "w");
			for (int i403 = 1; i403 <= nb_task - i; ++i403)
				fprintf(f, "%d %d %d\n", i403, ed[node->tache[i + i403 - 1]] - t0, pt[node->tache[i + i403 - 1]]);
			fclose(f);
			spawnl(P_WAIT, "chudc.exe", "chudc.exe", "d.dat", NULL);
			f = fopen("chudc.txt", "r");
			int op = -1; double tmp;
			fscanf(f, "%lf\n", &tmp);
			fscanf(f, "%d\n", &op);
			fclose(f);
			int trueSol = op + (nb_task - i)*t0;
			//if (isLb && sol > trueSol || !isLb && sol != trueSol){
			if ( !isLb && sol != trueSol){
				printf("isLB=%d; %d!=%d. t0=%d; v->niv=%d; v->temps=%d;\n",isLb, sol, trueSol, t0, i, node->temp);
				getch();
				sol = trueSol;
			}
		}
		//if (!isLb && t0 == 480 && i == 3)			printf("hehe");

		if (FALSE &&!isLb ){
			static int inspt[] = { 41, 62, 52, 53, 0, 22, 39, 29, 42, 40, 48, 57, 35, 2, 34, 23, 65, 68, 8, 58, 26, 63, 45, 3, 44, 54, 13, 14, 12, 5, 25, 37, 15, 32, 16, 38, 46, 49, 27, 17, 24, 28, 18, 50, 6, 10, 21, 7, 9, 43, 30, 33, 55, 1, 4, 31, 20, 36,
				 11, 19, 36, 51, 56, 59, 60, 61, 64, 66, 67, 69 };
			static int m = 0, i414 = 0;
			PTR*v = node;
			int flag[100] = { 0 };
			for (i414 = i; v != NULL && i414 < nb_task ; ++i414)flag[v->tache[i414]] = 1;
			for (i414 = 69; v != NULL && i414 >=0; --i414){
				if (flag[inspt[i414]] != 1)
				//if (v->tache[i414] != inspt[i414])
					break;
			}
			if ((69-i414) == nb_task-i)
			//if (i414 + 1 > m)
			//{
			//	m = i414 + 1;
			//	printf("%d\n", m);
			//}
			//if (i414 == v->niv)
			{
				//maxniv = i414;
				//if (maxniv == 58)
				printf("%d\t%d\t%d\t%d\n", 57-i414, t0, sol, isLb);
				FILE*f = fopen("d.dat", "w");
				for (int i403 = 1; i403 <= nb_task - i; ++i403)
					fprintf(f, "%d %d %d\n", i403, ed[node->tache[i + i403 - 1]] - t0, pt[node->tache[i + i403 - 1]]);
				fclose(f);
				spawnl(P_WAIT, "chudc.exe", "chudc.exe", "d.dat", NULL);
				f = fopen("chudc.txt", "r");
				int op = -1; double tmp;
				fscanf(f, "%lf\n", &tmp);
				fscanf(f, "%d\n", &op);
				fclose(f);
				int trueSol = op + (nb_task - i)*t0;
				//if (isLb && sol > trueSol || !isLb && sol != trueSol){
				if (!isLb && sol != trueSol){
					printf("isLB=%d; %d!=%d. t0=%d; v->niv=%d; v->temps=%d;\n", isLb, sol, trueSol, t0, i, node->temp);
					getch();
					sol = trueSol;
				}
			}
		}
		//if (isLb)printf("add\n");
		DBAddPb(node->tache + i, nb_task - i , t0, sol, isLb);	//! Note we do not store sequence. Only store jobset
	}
	if (sol + Ci < bestlb[node->niv - 1])	
		bestlb[node->niv - 1] = sol + Ci; //! update father's bestlb: opt is also lb
	if (!isLb && sol + Ci < best[node->niv - 1])	best[node->niv - 1] = sol + Ci; //! update father's best
	if (memlb[node->niv] == 1) memlb[node->niv - 1] = 1;//!raise the pollution
	best[node->niv] = INT_MAX;					//! Empty levels should have their best value reset
	bestlb[node->niv] = INT_MAX;
	memlb[node->niv] = 0;
}
double toptcpu;
//long long maxram=0;
long long get_ram_usage();
double get_cpu_time();
double get_wall_time();
//========== EoS ==========


// SearchStrategy:
// 0, classical depth first
// 1, classical best first
// 2, classical breadth first
// 3, depth first memo
// 4, best first memo
// 5, breadth first memo
// 6, solution/lb memo (depth first)
int SearchStrategy = -1;	//! put -1 to read from ini

//! DB cleaning strategy: (just use 3)
// 0, default=shortest first; 
// 1, LUFO_old (new item always has its counter = 0); 
// 2, double longest first, multiple deletion;
// 3, LUFO (new item's counter value considers the counter of nodes that it dominated)
int CleanStrategy = 3;

//! Dominance Strategy
// 12, enable predictive node cutting.
// 0, disable
// Details: The last 4 bits: whether add k perm found to db;  GenerateKPerm; Whether Generate all kperm and add the best (otherwise add the first one); Enable Dom cond on unscheduled jobs (Jouglet);
int DomStrategy = 0; 
extern int TimesClean;
extern long long NbCleanMin, NbCleanAvg, NbCleanMax;
extern int CutActive, CutDone;
extern int CutBef;
extern long long NbAddKPerm;
optb(s)
char *s;
{
	int i, br;
	double x1;
	double x2;
	FILE *fichier;
	float tpsbsr;

	for (i = 0; i<nb_task; ++i) premier[i] = NULL;
	//lecture("85.37.23.22.21.20.dat");
	//lecture("d.dat");
	//lecture("85.txt");
	//lecture("10/1.txt");
	//lecture("70.66.txt");
	//lecture("../../Output/data/70/20.txt");
	lecture(s);
	bsup = INFINI;
	for (i = 1; i<5; ++i){
		br = sch(i, 0, nb_task, buf->tache);
		if (br<bsup)	bsup = rsc = br;
	}
	//SearchStrategy = 6; printf("Stra=0 fixed !!!\n");
	if (SearchStrategy>2)	AllocDB(DBMaximumDimension, nb_task);
	if (SearchStrategy == 6){
		DomStrategy = 1;
		printf("DomStrategy set to 1 for Strategy 6. Dom based on unscheduled jobs\n");
		for (i = 0; i < NB + 1; ++i){
			best[i] = INT_MAX;
			bestlb[i] = INT_MAX;
			memlb[i] = 0;
		}
		//printf("70/6.txt should = 88787\n");
		//printf("70/12.txt should = 81355\n");
		//printf("70/19.txt should = 78384\n");
		//printf("70.35txt should = 104710\n");
		//printf("70.66.txt should = 104144\n"); 
		//printf("70.128.txt should = 130889\n");
		//printf("70.149.txt should = 167249\n");
		printf("70.38.txt should = 90227\n");
	}
	bsup -= somrel;
	//maxram = get_ram_usage();
	x1 = get_wall_time();
	x2 = get_cpu_time();

	if (nb_task <= NB_MAX) solopt(s);
	//rsc-=somrel;
	topt = (get_wall_time() - x1) ;
	toptcpu = get_cpu_time() - x2;
}

lecture(f_name)
char *f_name;
{
	printf("Reading %s.\n", f_name);
	FILE *f;
	int i, ordjob[NB], num;	//! ordjob = job id, start from 0
	Echelle hie;
	PTR *affptr();
	nb_task = 0;
	if ((f = fopen(f_name, "r")) == NULL)
	{
		printf("There are problems to open the file %s\n", f_name);
		exit(0);
	}
	else
	{
		for (somrel = 0, i = 0; fscanf(f, "%d %d %d\n",
			&num, &ed[i], &pt[i]) != EOF; ++i)
		{
			somrel += ed[i];
			nb_task++;
		}
		/*for(somrel=0,i=0;i<nb_task && fscanf(f,"%d %d %d\n",
		&num,&ed[i],&pt[i])!=EOF;++i) somrel+=ed[i];*/
		fclose(f);
	}
	buf = affptr();
	buf->niv = buf->brn = 0;
	buf->temp = 0;
	hie.ncrt = 2;
	hie.signe = (short *)malloc(2 * sizeof(short));
	hie.tend = (short *)malloc(2 * sizeof(short));
	hie.crit = (PT *)malloc(2 * sizeof(PT));
	for (i = 0; i<2; ++i)
	{
		*(hie.signe + i) = vrai;	//! int
		*(hie.tend + i) = faux;		//! increasing
	}
	hie.crit->pint = pt;
	(hie.crit + 1)->pint = ed;
	trie(faux, (int)nb_task, (int *)NULL, hie);//! sort hie to spt then on ert(earliest release date)
	//!LOG
	LOG("Instance read. N : %d\nJobs : \n", nb_task);
	LOG_FOR(j, nb_task, "%d\t", j); LOG("\n");
	LOG_FOR(j, nb_task, "%d\t", ed[j]); LOG("\n");
	LOG_FOR(j, nb_task, "%d\t", pt[j]);
	LOG("\nSpt : \n");
	LOG_FOR(j, nb_task, "%d\t", (hie.crit[0]).pint[j]);

	for (i = 0; i<nb_task; ++i) ordjob[i] = i;
	hie.crit->pint = ed;
	(hie.crit + 1)->pint = ordjob;
	trie(vrai, (int)nb_task, buf->tache, hie); //! sort hie to ert then spt. Hie is not modified, the order is put in buf->tache. So the content of tache is the order of ert

	//!LOG
	LOG("\nOrderErt : \n");
	LOG_FOR(j, nb_task, "%d\t", buf->tache[j]);

	free((char *)hie.crit);
	free((char *)hie.signe);
	free((char *)hie.tend);

	for (i = 0; i<NBBIT; ++i) buf->stat[i] = 0;
	// We read the search strategy //! and the db cleaning strategy
	//!
	if (SearchStrategy < 0){
		f = fopen("chudc.ini", "rt");
		fscanf(f, "%d\n", &SearchStrategy);
		fscanf(f, "%d\n", &CleanStrategy);
		fscanf(f, "%d\n", &DomStrategy);
		fclose(f);
	}
	// If SearchStrategy=0 then we have a classical depth first
	// If SearchStrategy=1 then we have a classical best first
	// If SearchStrategy=2 then we have a classical breadth first
	// If SearchStrategy=3 then we have a classical depth first + dc
	// If SearchStrategy=4 then we have a classical best first + dc
	// If SearchStrategy=5 then we have a classical breadth first + dc
	printf("\nUsing file '%s'; %d jobs; strategy %d, cleanstra %d, domstrategy %d.\n", f_name, nb_task, SearchStrategy, CleanStrategy, DomStrategy);
}

int nodecompare(PTR *first, PTR *second)
{ // Returns 0 if firts dominates second
	// Returns 1 if none dominates
	int r, i;

	// We first test using the condition of Chu
	r = nb_task - second->niv;
	if (first->brn <= second->brn && r*R(first->temp, first->tache[first->niv]) + first->brn
		<= r*R(second->temp, second->tache[second->niv]) + second->brn)
	{
		for (i = second->niv; i<nb_task && second->tache[i] == first->tache[i]; ++i);
		if (i == nb_task) return(0);
	}

	// We next test using the condition of Chand et al.
	/*if ((second->brn>first->brn)&&(((signed int)second->brn-(signed int)first->brn)>=(signed int)(r-1)*((signed int)pt[second->tache[second->niv-1]]-(signed int)pt[first->tache[first->niv-1]])))
	{
	for(i=second->niv;i<nb_task && second->tache[i]==first->tache[i];++i);
	if(i==nb_task) return(0);
	}
	if ((second->brn<first->brn)&&(((signed int)second->brn-(signed int)first->brn)*(signed int)r<=((signed int)pt[second->tache[second->niv-1]]-(signed int)pt[first->tache[first->niv-1]])))
	{
	for(i=second->niv;i<nb_task && second->tache[i]==first->tache[i];++i);
	if(i==nb_task) return(0);
	}*/
	if ((second->brn>first->brn) && (((signed int)second->brn - (signed int)first->brn) >= (signed int)r*((signed int)first->temp - (signed int)second->temp)))
	{
		for (i = second->niv; i<nb_task && second->tache[i] == first->tache[i]; ++i);
		if (i == nb_task) return(0);
	}
	return(1);
}

//! Return 1 if j dominates k. On unscheduled jobs, not polluting.
int domJouglet(int k, int j, int maxrNS, int niv){
	int t1, t2, t3, rk = ed[k], pk = pt[k], rj = ed[j], pj = pt[j];
	//! by cases
	if (rj + pj <= rk + pk && pj < pk && maxrNS >= rk + pk){
		//case1
		t1 = 0; t2 = t3 = pk - pj;
	}
	else if (rj + pj <= rk + pk && pj < pk && maxrNS < rk + pk  && rj <= rk && rk + pj < maxrNS){
		//case2
		t1 = maxrNS - rk - pk;
		t2 = t3 = maxrNS - rk - pj;
	}
	else if (rj + pj <= rk + pk && pj < pk&& maxrNS < rk + pk && rj <= rk &&maxrNS <= rk + pj){
		//case3
		t1 = max(rj + pj, maxrNS) - rk - pk; t2 = t3 = 0;
	}
	else if (rj + pj <= rk + pk && pj < pk && maxrNS < rk + pk && rj>rk){
		//case4
		t1 = max(rj + pj, maxrNS) - rk - pk;
		t2 = t3 = max(rj + pj, maxrNS) - rk - pj;
	}
	else if (rj + pj <= rk + pk && pj >= pk && maxrNS >= rk + pk){
		//case5
		t1 = t2 = t3 = 0;
	}
	else if (rj + pj <= rk + pk && pj >= pk && maxrNS < rk + pk){
		//case6
		t1 = max(rj + pj, maxrNS) - rk - pk; t2 = t3 = 0;
	}
	else if (rj + pj > rk + pk && pj < pk){
		//case7
		t1 = rj + pj - rk - pk; t2 = t3 = rj - rk;
	}
	else if (rj + pj > rk + pk && pj > pk){
		//case8
		t1 = t2 = rj + pj - rk - pk; t3 = max(0, rj - rk);
	}
	//! Gamma(t,t1,t2,t3)
	int gamma, t = rk + pk;
	int wj = 1, wk = 1, dj = 0, dk = 0; //! do not modify
	//! for gamma'
	gamma = wj*max(0, t + pj - dj) - wk*max(0, t + t1 + pk - dk)
		+ wk*max(0, rk + pk - dk) - wj*max(0, rj + pj - dj)
		- t2* (nb_task - niv - 2) - 0;
	if (gamma > 0) return 1;	//! I avoided put >=

	//! for Gamma
	t = rk + pk;
	t1 = rj + pj - rk - pk;
	t2 = rj + pj - rk;
	t3 = max(0, rj - rk);
	gamma = wj*max(0, t + pj - dj) - wk*max(0, t + t1 + pk - dk)
		+ wk*max(0, rk + pk - dk) - wj*max(0, rj + pj - dj)
		- t2* (nb_task - niv - 2) - 0;

	if (gamma > 0) return 1;    //! I avoided put >=
	return 0;
}

solopt(char* datafile)
{
	int i, j, nb, bm, t, t1, idle, w, w1, w2, n, C, mu, x, *med, cut, index2;
	PTR *v, *parc, *tmp, *pparc, *preclast;

	nbatt1 = nbatt = nbpas = nbsepa = nbPolluted= 0;
	//printf("This is a tmp version for testing.\n Use srptimp for LB. Set ub=518736. [Optimp.c ligne 426.]\n");
	srpt(buf); //should be enabled in the normal version
	//srptImp(buf, datafile); //! just for testing
	//rsc = 518736;			//! just for testing

	LBroot = buf->binf;
	
	/*
	int ub = rsc;
	int nbStep = 4;
	if (rsc - LBroot < 4) nbStep = 1;
	double stepUB = (rsc - LBroot) / (double)(nbStep);
	for(int iStepUB=1, rsc = LBroot + stepUB; iStepUB <= nbStep; rsc=LBroot+stepUB*iStep)

	//! tentative ub
	//for each iteration
	for (i = 0; i<nb_task; ++i) premier[i] = NULL;
	last == NULL;
	//while (preclast->next != last) preclast = preclast->next;
	//preclast = premier[0];
	//while (preclast->next != last) preclast = preclast->next;*/

	LOG_NODE(buf);
	//LOG("\nrsc=%d\n", rsc);
	if (buf->binf<rsc)	sauve(buf);
	if (buf->binf == rsc) ++nbsepa;
	for (;;)
	{
		// We extract a node
		if ((SearchStrategy == 0) || (SearchStrategy == 3) || (SearchStrategy == 6) || SearchStrategy==7)
		{ // Depth first strategy. //!Last is the last node of the whole node list of the tree
			if (last == NULL)
			{
				preclast = NULL;
				last = premier[0];
			}
			if (last != NULL)
			{
				if (last->next != NULL)
				{
					while (last->next != NULL)
					{
						preclast = last;
						last = last->next;
					}
				}
				else
				{ // We have to recompute preclast
					if (last == premier[0]) preclast = NULL;
					else
					{
						preclast = premier[0];
						while (preclast->next != last) preclast = preclast->next;
					}
				}
				nb = 0;
			}
			v = last;
			LOG("\nExtracted node : \n");
			LOG_NODE(v);
			if (preclast != NULL) preclast->next = last->next; //! that the node out of the list
			else premier[0] = NULL;
			last = preclast;
			cut = 0;
			if (v == NULL) nb = -1;

			if ((nb != -1) && (SearchStrategy == 3))
			{ // We compare with done nodes
				index2 = 100;
				if (v->niv>0 && DBAdd(v->tache, v->brn, v->temp, v->niv, 1, &index2) == -1)
				{ // The node v is dominated
					cut = 1;
					CutDone++;
				}

				if ((DomStrategy & 4) && cut == 0 && v->niv >= 3)
				{ // The current node is not dominated so we try to generate alternative sequences to try 
					// to find one dominating the current node : in this case it is added to the database
					// Added by VTkindt on the 04/03/2017
					cut = DBGenerate(v->tache, v->brn, v->temp, v->niv);
				}

				// If the current node is not dominated then we add it to the done nodes
				if (cut == 0)
				{
					// We save the done node
					if (v->niv>0)
						DBAdd(v->tache, v->brn, v->temp, v->niv, -1, NULL);
				}
			}
			//!MEMO is added here
			if ((nb != -1) && (SearchStrategy == 6))
			{ // We search in solved nodes

				if (v->binf < 0){//! This is a sentinel node. it is just solved and should be saved to database (solution memo)
					if (best[v->niv] != INT_MAX)//INTMAX when the node has no solution: all children are cut. we wanted to also memo this in order to indicate the node as not promising, however it is wrong since what is unpromising is not the subpb but the partial seq+subpb.
					{
						UpdateBestOnSolve(v, TRUE, FALSE);
					}
					else{
						if (bestlb[v->niv] == INT_MAX) //! if lb is not updated by children, remember the computed one.
							bestlb[v->niv] = -v->binf-1 - v->brn;
						UpdateBestOnSolve(v, TRUE, TRUE);	//! if some updates are missed, nodes may be cut innocently since the lb is higher than what it should be
					}
					cut = 1; //since this is a sentinel node
				}
				else if(v->niv>1){
					//!Find sol from db
					int sol = 999999;
					int t0 = 999999;
					t0 = ed[v->tache[v->niv]];
					if (t0 < v->temp) t0 = v->temp;
					int isLB = 0;
					if ((sol = DBSearchPb(v->tache + v->niv, nb_task - v->niv, t0, &isLB)) != -1){
						//! opt sol is also lb
						if (v->binf > sol + v->brn)
							printf("here!");
						v->binf = sol + v->brn;
						//assert(bestlb[v->niv] > sol); 
						bestlb[v->niv] = sol;
						//UpdateBestOnSolve(v, FALSE, TRUE);

						if (!isLB){
							if (sol + v->brn < rsc) {
								best[v->niv] = sol;
								rsc = sol + v->brn;
								//!also need to update its father. Only do this when sol<rsc since otherwise this is a cut node, which is not optimal. !!! Caution on this
								UpdateBestOnSolve(v, FALSE, FALSE);
							}
							else{//not improving best, but it's still valid as a lb.
								//bestlb[v->niv] = sol;
								UpdateBestOnSolve(v, FALSE, TRUE);
							}
							cut = 1;
							CutDone++;
							LOG("\nSolMemo cut size %d\n", nb_task - v->niv);
						}
						else if (Stra6LBMemoOn){
							UpdateBestOnSolve(v, FALSE, TRUE);
							if (sol + v->brn > rsc){//Cut by UB
								cut = 1;
								CutActive++;
								LOG("\nLbMemo cut size %d\n", nb_task - v->niv);
							}
						}
					}
					else if (v->binf>=rsc){//cut by ub, which may have been updated after the creation of node v
						bestlb[v->niv] = v->binf - v->brn;
						UpdateBestOnSolve(v, TRUE, TRUE);
					}
				}
			}
		}
		if ((SearchStrategy == 1) || (SearchStrategy == 4))
		{ // Best first strategy
			for (bm = INFINI, nb = -1, i = 0; i<nb_task - 2; ++i)
				if (premier[i] != NULL && premier[i]->binf<bm)
				{
					nb = i;
					bm = premier[i]->binf;
				}
			if (nb != -1)
			{
				v = premier[nb];
				premier[nb] = premier[nb]->next;
			}
			cut = 0;
			if ((nb != -1) && (SearchStrategy == 4))
			{ // We compare with other nodes
				// We compare with active nodes
				pparc = v;
				parc = v->next;
				while ((cut == 0) && (parc != NULL))
				{
					if (nodecompare(parc, v) == 0)
					{
						cut = 1;
						CutActive++;
						//if (v->binf > parc->binf){
						//	printf("\nCurrent node:\n");
						//	logNode(v);
						//	printf("\nActive node:\n");
						//	logNode(parc);
						//	getchar();
						//}
						//Lei: no instance is found finally. But some cut are done between nodes in equality
					}
					else if (nodecompare(v, parc) == 0)
					{ // We remove parc
						tmp = parc;
						CutActive++;
						parc = parc->next;
						pparc->next = parc;
						if (premier[nb] == tmp) premier[nb] = parc;
						free((char *)(tmp->tache));
						free((char *)tmp);
					}
					else
					{
						pparc = parc;
						parc = parc->next;
					}
				}
				// We compare with done nodes		  	

				index2 = 100;
				if ((cut == 0) && v->niv>0 && (DBAdd(v->tache, v->brn, v->temp, v->niv, 1, &index2) == -1))
				{ // The node v is dominated
					cut = 1;
					CutDone++;
				}
				if ((DomStrategy & 4) && cut == 0 && v->niv >= 3)
				{ // The current node is not dominated so we try to generate alternative sequences to try 
					// to find one dominating the current node : in this case it is added to the database
					// Added by VTkindt on the 04/03/2017
					cut = DBGenerate(v->tache, v->brn, v->temp, v->niv);
				}
				// The node is not dominated so we add it to the list of done nodes
				if (cut == 0)
				{
					// We save the done node
					if (v->niv>0)
						DBAdd(v->tache, v->brn, v->temp, v->niv, -1, NULL);
				}
			}
		}
		if ((SearchStrategy == 2) || (SearchStrategy == 5))
		{ // Breadth first strategy
			for (bm = INFINI, nb = -1, i = 0; i<nb_task - 2; ++i)
				if (premier[i] != NULL)
				{
					nb = i;
					break;
				}
			cut = 0;
			if (nb != -1)
			{
				v = premier[nb];
				premier[nb] = premier[nb]->next;
			}
			if ((nb != -1) && (SearchStrategy == 5))
			{ // We compare with active nodes
				pparc = v;
				parc = v->next;
				cut = 0;
				while ((cut == 0) && (parc != NULL))
				{
					if (nodecompare(parc, v) == 0)
					{
						cut = 1;
						CutActive++;
					}
					else if (nodecompare(v, parc) == 0)
					{ // We remove parc
						tmp = parc;
						CutActive++;
						parc = parc->next;
						pparc->next = parc;
						if (premier[nb] == tmp) premier[nb] = parc;
						free((char *)(tmp->tache));
						free((char *)tmp);
					}
					else
					{
						pparc = parc;
						parc = parc->next;
					}
				}
			}
		}
		if (nb == -1) break;
		if (cut == 0)
		{ // The current node has not been cut
			LOG("\nNode not cut.");
			//if (v->niv==1 && v->brn==219)printf("not ok");
			if (v->binf<rsc){
				BOOL isSentinelFatherSaved = FALSE;
				for (n = v->niv, i = v->niv + 1; i<nb_task; ++i) //! tache contains spt job id, the order in taches is in ert
					if (v->tache[i]<v->tache[n]) n = i;
				x = v->tache[n];			//! x = job with spt
				C = C1[n] = F(v->temp, x);  //! finishing time of the spt job 
				etat[n] = 1;
				mu = nb_task - v->niv;		//! mu=Nb jobs to scheduling
				for (i = v->niv; R(v->temp, (w = v->tache[i])) + pt[x] < C; ++i)
				{
					if ((C1[i] = F(v->temp, w)) < C) C = C1[i];	//! C, earliest finishing time
					etat[i] = 1;
					//if (SearchStrategy != 6)
					{//! dominance conditions are disabled for stra6
						for (t = 0, bm = v->niv + 1, j = 0; j < v->niv; ++j)//! using the fixed jobs, therefore polluting
						{
							idle = F(t, v->tache[j]);
							t1 = pt[v->tache[j]] - pt[w];
							if (((w1 = idle - F(t, w)) > 0 || w1 == 0 &&
								v->tache[j] > w) && w1 >= mu*t1  /*Condition 9 du papier de chu*/
								||
								w1<0 && bm*w1 >= t1 /* Condition 10 du papier de Chu*/
								|| R(t, v->tache[j]) >= R(t, w) &&
								v->temp - (bm - 1)*(F(t, w) - R(t, v->tache[j]))>R(t, w)) /* ???? */
							{
								etat[i] = 0;	//! eliminating job candidates
								if (SearchStrategy == 6) etat[i] = -1;//! should evaluate and then cut
								break;
							}
							t = idle;
							--bm;
						}
					}
				}
				if (C<C1[n]) n = i - 1;
				for (x = i; i<n; ++i) etat[i] = faux; //TODO could here be a dom condition?No, only when it involves fixed jobs
				int lastCandidate = n;
				//! find max r among unscheduled jobs
				//				for (; lastCandidate >= v->niv; --lastCandidate)
				//				if (etat[lastCandidate]) break;
				for (j = v->niv; j <= lastCandidate; ++j){
					//! Lei: Add jouglet's conditions here: 13.8.3 based on unscheduled jobs
					if (etat[j] && (DomStrategy & 1)){
						int k = j, jj; //! To fit the var names in the theorems
						for (jj = k + 1; jj <= lastCandidate; ++jj)
							if (etat[jj]==1 && domJouglet(v->tache[k], v->tache[jj], ed[v->tache[lastCandidate]], v->niv)){
								etat[k] = 0;
								//if (SearchStrategy == 6) etat[k] = -1; //! which means "evaluate then cut"
								//!TODO whether add the dominant one to the db ? No, since the dominant one will be branched on and added to db automatically
								break;
							}
						//! End of jouglet's conditions
					}
					if (etat[j])	//! branching on candidate jobs. could =-1
					{
						//if (SearchStrategy != 6){//! seems not based on fixed jobs, can be kept for stra6
							w = pt[v->tache[j]];
							for (i = v->niv; i < j && i < x; ++i){
								t = w - pt[v->tache[i]];
								if ((w2 = C1[j] - C1[i]) >= 0 && w2 >= t*(mu - 1) /* Theorem 6 of Chu */
									|| w2 <= 0 && w2*mu >= t)  /* Theorem 7 of chu */
									break;
							}
						//}
						//etat[j] = (i == j || i == x);
						//if (etat[j] == 0 && SearchStrategy == 6)etat[j]=
						if (i == j || i == x)//! disable dom for stra6
						{
							med = buf->tache;
							*buf = *v;
							buf->tache = med;
							memcpy((char *)(buf->tache),
								(char *)(v->tache), nb_task*sizeof(int));
							if (buf->binf < 0)//! since we modified v->inf in this loop......血与泪
								buf->binf = -buf->binf;
							if (isSentinelFatherSaved == FALSE && (SearchStrategy==6 || SearchStrategy==7)){
								isSentinelFatherSaved = TRUE;
								//if(SearchStrategy==6) 
								best[v->niv] = INT_MAX; //! do this in stra7 will overwrite the lb value, computed when branching
								best[v->niv + 1] = INT_MAX;
								bestlb[v->niv] = INT_MAX;//v->binf-v->brn; //! 血与泪
								bestlb[v->niv + 1] = INT_MAX;
								v->binf = -v->binf-1; //! indicate this node as sentinel (explored father node)
								//!cancel effect of sauve
								--nbatt;
								--nbpas;
								--nbsepa;
								sauve(v);
							}
							separe(buf, j, etat[j]);
						}
					}
				}
			}
			//else { printf("Cut by ub, not expected.\n"); }
		}
		--nbatt;
		free((char *)(v->tache));
		free((char *)v);
	}// for true
}

separe(cop, x, etatx)
PTR *cop;
int x;
short etatx;
{
	int i, w, zi, wei;
	w = (*(cop->tache + x)); //! w is the job id to branch
	LOG("\nBranching on tache %d\n", w);
	//LOG_NODE(cop);
	for (i = x; i>cop->niv; --i) *(cop->tache + i) = (*(cop->tache + i - 1));
	*(cop->tache + cop->niv) = w;
	cop->temp = F(cop->temp, w);
	zi = w / NOI;
	wei = w%NOI;
	cop->stat[zi] += EXP2[wei];
	++cop->niv;
	cop->brn += cop->temp;							//! New LB for current node

	//if ( cop->tache[0] == 41 && cop->tache[1] == 62)		printf("here");
	if (etatx == -1 && cop->niv > 0) {
		memlb[cop->niv - 1] = 1;//! mark this node's father node as polluted. This node is not, since it might be solved directly to optimality.
		++nbPolluted;
	}
	if (FALSE && cop->niv == 58 && cop->temp==3444){
		static int inspt[] = { 41, 62, 52, 53, 0, 22, 39, 29, 42, 40, 48, 57, 35, 2, 34, 23, 65, 68, 8, 58, 26, 63, 45, 3, 44, 54, 13, 14, 12, 5, 25, 37, 15, 32, 16, 38, 46, 49, 27, 17, 24, 28, 18, 50, 6, 10, 21, 7, 9, 43, 30, 33, 55, 1, 4, 31, 20, 47 };
		static int m = 0, i414 = 0;
		PTR*v = cop;
		//int flag[100] = { 0 };
		//for (i414 = 0; v != NULL && i414 < v->niv && v->niv <= 58; ++i414)flag[inspt[i414]] = 1;
		for (i414 = 0; v != NULL && i414 < v->niv && v->niv <= 58; ++i414){
			//if (flag[v->tache[i414]] != 1)
			if (v->tache[i414] != inspt[i414])
				break;
		}
		//if (i414 == v->niv)
		////if (i414+1>m)
		//{
		//	//m = i414+1;
		//	printf("%d\n", m);
		//}
	}


	if (ed[cop->tache[nb_task - 1]] <= cop->temp){	//! Remaining jobs are all released, opt order = spt
		int oldBn = cop->brn;
		cop->brn += tspt(cop);
		if ((cop->brn) < rsc)//!Only update opt when rsc is improved, otherwise it's possible that the opt node was cut by dom condition and therefore not considered
		{
			rsc = cop->brn;
			
			//if (rsc == 167249){ printf("11111"); getch(); }
			if (etatx != -1){
				++nbsepa;
				LOG("1:");
				LOG_NODE(cop);
			}
			if (SearchStrategy == 6 || SearchStrategy==7) {
				bestlb[cop->niv] = best[cop->niv] = cop->brn - oldBn;
				UpdateBestOnSolve(cop, TRUE, FALSE);
			}
		}
		else if (SearchStrategy == 6 || SearchStrategy==7){
			// Not improving ub, still a good lb
			bestlb[cop->niv] = cop->brn - oldBn;
			UpdateBestOnSolve(cop, TRUE, TRUE);
		}
		LOG("\nBranched node solved directly.Return.");
		return;
	}
	if (cop->niv >= nb_task - 2)
	{//! The leave nodes are not saved: only ub are updated
		int oldBn = cop->brn;
		cop->brn += prtf(cop->temp, nb_task - cop->niv, cop->tache + cop->niv);
		if ((cop->brn) < rsc)//! add =, otherwise updateBest will rarely be called
		{
			rsc = cop->brn;
			//if (rsc == 167249){ 
				//printf("2222"); getch(); }
			if (etatx != -1){
				++nbsepa;
				LOG("2:");
				LOG_NODE(cop);
			}
			if (SearchStrategy == 6 || SearchStrategy == 7) {
				bestlb[cop->niv] = best[cop->niv] = cop->brn - oldBn;
				UpdateBestOnSolve(cop, TRUE, FALSE);
			}
		}
		else if (SearchStrategy == 6 || SearchStrategy == 7){
			// Not improving ub, still a good lb
			bestlb[cop->niv] = cop->brn - oldBn;
			UpdateBestOnSolve(cop, TRUE, TRUE);
		}
		LOG("\nBranched node is solved as leaf. Return.");
	}
	else{
		//! Here when trying to cut a node by LB, we do in a new way: search LB in DB, if not found then compute it and update best. When arrive at a sentinel node, save the best to db.
		//! BUGFIX: we can not do this in this code since the srpt() function also updates rsc, so it can not be skipped.
		//int sol = 999999;
		//int t0 = 999999;
		//int isLb=0;
		//if (SearchStrategy == 6 && cop->niv > 1 && Stra6LBMemoOn){
		//		//!Find sol from db
		//		t0 = ed[cop->tache[cop->niv]];
		//		if (t0 < cop->temp) t0 = cop->temp;
		//		if ((sol = DBSearchPb(cop->tache + cop->niv, nb_task - cop->niv, t0, &isLb)) != -1){
		//				cop->binf = sol + cop->brn;	//! we donot propagate this now: it will be done when the node is extracted.
		//		}
		//}
		//if (sol == 999999 || sol == -1)
		int updatedInSprt = FALSE;
		{
			int r = rsc;
			srpt(cop);
			//!!Bug fixed... srpt is also updating rsc..........So best/bestlb should be updated
			if (rsc < r){
				updatedInSprt = TRUE;
				bestlb[cop->niv] = best[cop->niv] = cop->binf - cop->brn;
				UpdateBestOnSolve(cop, TRUE, FALSE);
			}
			else if(etatx==-1){//! This node is to be cut, so the lb should be taken into account
				bestlb[cop->niv] = cop->binf - cop->brn;
				UpdateBestOnSolve(cop, TRUE, TRUE);
			}
		}
		if (cop->binf<rsc)
		{
			if ((zi = cop->brn + aprtf(cop->temp, nb_task - cop->niv,
				cop->tache + cop->niv))<rsc)
				rsc = zi;
			if (etatx == 1){
				sauve(cop);
				LOG("3:");
				LOG_NODE(cop);
			}
			LOG("\nBranched node is saved to the tree.");
		}
		else if(SearchStrategy==6){
			// The node is cut by lb memo. We may want to update its father
			//if (sol < 999999 && sol > -1)
			//{
			//	CutActive++;
			//	/*if (!isLb ){ //! bug fixed: should not have updated best since rsc is not improved!
			//		best[cop->niv] = bestlb[cop->niv] = sol;
			//		UpdateBestOnSolve(cop, FALSE, FALSE);
			//	}
			//	else*/
			//	{
			//		bestlb[cop->niv] = sol;
			//		UpdateBestOnSolve(cop, FALSE, TRUE);
			//	}
			//} else 
			if (!updatedInSprt)
			{
				// We can add this binf since it has no child. //! bug fixed: for lb, no end nodes should be omitted for the update of lb. Otherwise lb may be greater than it should and the node be cut by fault
				bestlb[cop->niv] = cop->binf - cop->brn;
				UpdateBestOnSolve(cop, TRUE, TRUE);
			}
			LOG("\nBranched node is cut by UB.");
		}
	}
}

tspt(cop)
PTR *cop;
{
	int i, t, ret;

	for (ret = 0, t = cop->temp, i = 0; i<nb_task; ++i){
		if ((cop->stat[i / NOI] / EXP2[i%NOI]) % 2 == 0){
			t += pt[i];
			ret += t;
		}
	}
	return(ret);
}
//! Heuristic to compute LB of the whole node
srpt(cop)
PTR *cop;
{
	int i, n, t, x, lq, cod, aj, t1;

	n = nb_task - cop->niv;
	lq = 0;
	for (i = 0; i<n; ++i){
		arb[i] = -1;
		qu[i] = 0;
		arb1[i] = 1000000000;
	}
	for (i = 0; i<n; ++i){
		deb1[i] = R(cop->temp, cop->tache[i + cop->niv]);
		prov1[i] = pt[cop->tache[i + cop->niv]];
	}
	for (aj = 0, cod = vrai, t = cop->temp, cop->binf = cop->brn, x = -1, i = 0; i<n;){
		if (deb1[i] + prov1[i]<t)
		{
			cod = faux;
			prov1[x] = t - deb1[i];
			t = deb1[i] + prov1[i];
			ajarb(i, prov1, &lq, arb, qu);
			x = i++;
		}
		else if (deb1[i] <= t)
		{
			if (x != -1) prov1[x] = t - deb1[i];
			ajarb(i, prov1, &lq, arb, qu);
			++i;
		}
		else {
			if (x != -1)
			{
				modarb(n, prov1, &lq, arb, qu);
				cop->binf += t;
			}
			if ((x = arb[0]) == -1)	t = deb1[i];
			else	
				t += prov1[x];
		}
	}
	if (x != -1){
		cop->binf += t;
		modarb(n, prov1, &lq, arb, qu);
	}
	for (; (x = arb[0]) != -1; modarb(n, prov1, &lq, arb, qu))
	{
		t += prov1[x];
		cop->binf += t;
	}
	if (cop->niv == 0 && cop->binf>rsc)
		printf("Erreur\n");
	if (cod && cop->binf<rsc)
		rsc = cop->binf;
}

//! Call the exe to get a good lb using srpt improved.
srptImp(PTR *cop, const char* datafile)
{
	printf("Using improved SRPT. Input=%s\n", datafile);
	float temps;
	system("del lbsrptimp.txt");
	//_spawnl(P_WAIT, "lbsrptimp.exe", "lbsrptimp.exe", datafile, NULL);
	system("lbsrptimp.exe donnees.dat");
	FILE * fichier = fopen("lbsrptimp.txt", "r");
	fscanf(fichier, "%lf\n", &temps);
	//printf("LB=%d\n", cop->binf);
	fscanf(fichier, "%d\n", &cop->binf);
	fclose(fichier);
	printf("LB=%d\n", cop->binf);
}

PTR *affptr()
{
	PTR *v;

	if ((v = (PTR *)malloc(sizeof(PTR))) == (PTR *)NULL) exiterr();
	if ((v->tache = (int *)malloc(nb_task*sizeof(int))) == (int *)NULL) exiterr();
	return(v);
}
//! Save node to the tree node list end. TO enable memo, it may be necessary to add new node to the end instead of a position related to its lb
sauve(cop)
PTR *cop;
{
	PTR *v, *v1, *vs;
	int i, *w, r, vn;

	if (SearchStrategy == 0 || SearchStrategy == 3 || SearchStrategy == 6 || SearchStrategy==7)
	{
		v1 = last;
		if (v1 == NULL) v = premier[0];
		else v = last->next;
		//if (SearchStrategy == 6 || SearchStrategy==7){
		//	for (; v != NULL //&& v->binf >= cop->binf //! for memo, we want the new node is exactly added at the end, since we are adding sentinel nodes (father nodes)
		//		; v = v->next) //!smaller LB, later in the list, earlier extracted
		//	{
		//		//printf("line 749 reached\n");
		//		v1 = v;
		//	}
		//}
		//else{
			for (; v != NULL && (v->binf >= cop->binf || v->binf<0); v = v->next) //!smaller LB, later in the list, earlier extracted. At the end inf increasing: v1->cop->v
			{
				v1 = v;
			}
		//}
	}
	else
	{
		v1 = NULL;
		v = premier[cop->niv];
		for (; v != NULL && v->binf <= cop->binf; v = v->next)
		{
			v1 = v;
		}
	}
	vs = affptr();
	w = vs->tache;
	*vs = *cop;
	vs->tache = w;
	memcpy((char *)(vs->tache), (char *)(cop->tache), nb_task*sizeof(int));
	vs->next = v;
	if (v1 == NULL)
	{
		if (SearchStrategy == 0 || SearchStrategy == 3 || SearchStrategy == 6 || SearchStrategy==7) premier[0] = vs;
		else premier[cop->niv] = vs;
	}
	else v1->next = vs;
	++nbatt;
	//printf("\rNbre de noeuds actifs : %ld",nbatt);
	if (nbatt>nbatt1) nbatt1 = nbatt;
	++nbpas;
	++nbsepa;
}

exiterr()
{
	printf("Probleme d'affectation Memoire, nbatt=%d\n", nbatt);

	exit(0);
}
//! Used to compute UB
aprtf(t, n, job)
int t, n, *job;
{
	int i, jw, ft, x, lq, a, b, f, c, c1, c2, e, t1;

	lq = 0;
	for (i = 0; i<n; ++i)
	{
		arb[i] = -1;
		qu[i] = 0;
	}
	for (i = 0; i<n; ++i)
	{
		deb1[i] = R(t, job[i]);
		prov1[i] = pt[job[i]];
	}
	for (a = n - 2, ft = 0, i = 0, x = -1; i<n;)
	{
		if ((b = arb[0]) == -1)
		{
			t1 = INFINI;
			b = i;
		}
		else{
			modarb(n, prov1, &lq, arb, qu);
			if (arb[0] != -1)	t1 = t;
			else	t1 = INFINI;
		}
		c2 = R(t, job[b]);
		f = c1 = 2 * c2 + prov1[b];
		c2 += prov1[b];
		for (jw = -1, i = max(b + 1, i); i<n && deb1[i]<c2; ++i)
		{
			if (deb1[i] <= t && prov1[i]<prov1[b])
			{
				ajarb(b, prov1, &lq, arb, qu);
				if (deb1[b]<t1)	t1 = deb1[b];
				c1 = f = f - prov1[b] + prov1[i];
				b = i;
			}
			if ((e = F(t, job[i]))<c2) c2 = e;
			if ((c = 2 * e - prov1[i])<c1)
			{
				if (jw != -1)
				{
					ajarb(jw, prov1, &lq, arb, qu);
					if (deb1[jw]<t1)	t1 = deb1[jw];
				}
				c1 = c;
				jw = i;
			}
			else	if (i != b)
			{
				ajarb(i, prov1, &lq, arb, qu);
				if (deb1[i]<t1)	t1 = deb1[i];
			}
		}
		if (jw != -1 && jw<n - 1 && deb1[jw + 1]<t1) t1 = deb1[jw + 1];
		if (jw == -1) x = b;
		else	if (f - c1 >= a*min(deb1[jw] - R(t, job[b]), deb1[jw] + prov1[jw] + prov1[b] - t1))
		{
			ajarb(b, prov1, &lq, arb, qu);
			x = jw;
		}
		else{
			x = b;
			ajarb(jw, prov1, &lq, arb, qu);
		}
		/*		if(n==30)
		printf("%d\t",R(t,job[x]));*/
		t = F(t, job[x]);
		--a;
		/*		if(n==30)
		printf("%d\t%d\n",t,job[x]+1);*/
		ft += t;
	}
	for (; (x = arb[0]) != -1; modarb(n, prov1, &lq, arb, qu))
	{
		/*		if(n==30)
		printf("%d\t%d\t%d\n",t,t+prov1[x],job[x]+1);*/
		t += prov1[x];
		ft += t;
	}
	return(ft);
}
//! Algo related to Priority Rule. See the paper
prtf(t, n, job)
int t, n, *job;
{
	int i, ft, x, lq, c, c1;

	lq = 0;
	for (i = 0; i<n; ++i)
	{
		arb[i] = -1;
		qu[i] = 0;
	}
	for (i = 0; i<n; ++i)
	{
		deb1[i] = R(t, job[i]);
		prov1[i] = pt[job[i]];
	}
	for (ft = 0, i = 0; i<n;)
	{
		if ((x = arb[0]) == -1) x = i;
		else	modarb(n, prov1, &lq, arb, qu);
		c1 = 2 * R(t, job[x]) + prov1[x];
		for (i = max(x + 1, i); i<n && (c = 2 * R(t, job[i]))<c1; ++i)
			if ((c = c + prov1[i])<c1)
			{
				ajarb(x, prov1, &lq, arb, qu);
				c1 = c;
				x = i;
			}
			else	ajarb(i, prov1, &lq, arb, qu);
			t = F(t, job[x]);
			ft += t;
	}
	for (; (x = arb[0]) != -1; modarb(n, prov1, &lq, arb, qu))
	{
		t += prov1[x];
		ft += t;
	}
	return(ft);
}

est(t, n, job)
int t, n, *job;
{
	int i, ft, x, lq;

	lq = 0;
	for (i = 0; i<n; ++i)
	{
		arb[i] = -1;
		qu[i] = 0;
	}
	for (i = 0; i<n; ++i) prov1[i] = pt[job[i]];
	for (ft = 0, i = 0; i<n;)
		if (ed[job[i]] <= t)
		{
			ajarb(i, prov1, &lq, arb, qu);
			++i;
		}
		else	if ((x = arb[0]) == -1) t = ed[job[i]];
		else
		{
			modarb(n, prov1, &lq, arb, qu);
			t += prov1[x];
			ft += t;
		}
		for (; (x = arb[0]) != -1; modarb(n, prov1, &lq, arb, qu))
		{
			t += prov1[x];
			ft += t;
		}
		return(ft);
}

ect(t, n, job)
int t, n, *job;
{
	int i, ft, x, lq, c;

	lq = 0;
	for (i = 0; i<n; ++i)
	{
		arb[i] = -1;
		qu[i] = 0;
	}
	for (i = 0; i<n; ++i)
	{
		deb1[i] = R(t, job[i]);
		prov1[i] = pt[job[i]];
	}
	for (ft = 0, i = 0; i<n;)
	{
		if ((x = arb[0]) == -1) x = i;
		else	modarb(n, prov1, &lq, arb, qu);
		t = F(t, job[x]);
		for (i = max(x + 1, i); i<n && deb1[i]<t; ++i)
			if ((c = deb1[i] + prov1[i])<t)
			{
				ajarb(x, prov1, &lq, arb, qu);
				t = c;
				x = i;
			}
			else	ajarb(i, prov1, &lq, arb, qu);
			ft += t;
	}
	for (; (x = arb[0]) != -1; modarb(n, prov1, &lq, arb, qu))
	{
		t += prov1[x];
		ft += t;
	}
	return(ft);
}

init(n, job)
int	n, *job;
{
	int	i;

	for (i = 0; i<n; ++i)
	{
		arb[i] = -1;
		qu[i] = 0;
	}
	for (i = 0; i<n; ++i)
	{
		deb1[i] = job[i];
		prov1[i] = pt[deb1[i]];
	}
	return(0);
}

sch(x, t, n, job)
int x, t, n, *job;
{
	switch (x)
	{
	case 1: return(aprtf(t, n, job));
	case 2: return(prtf(t, n, job));
	case 3: return(est(t, n, job));
	case 4:	return(ect(t, n, job));
	default:printf("Il y a des erreurs dans votre programme!!!\n");
		return(INFINI);
	}
}

//#define pt 0.00001

main(argc, argv)
int argc;
char *argv[];
{
	setvbuf(stdout, NULL, _IONBF, 0);
	int i;
	double nnode, nsepa, nnpas, temps, someh, somop;
	char numex[8];
	FILE *fichier;
	time_t now;
	time(&now);
	printf(" ===== In solver : %s", ctime(&now));
	printf("MaxJobs set to 220, DbDimention set to 550w.\n");
	/*if(argc==1)
	{
	printf("Erreur: pas d'arguments!!!\n");
	exit(0);
	}*/
	//printf("Exemple\tNb_Prd\tSepar\tNodgen\tNodes\tTime\tHeur\tOpt\n");
	for (nnpas = nnode = nsepa = temps = someh = somop = 0.0, i = 1; i<2; ++i)
	{
		//replch(numex,"donnees.dat");
		//strcpy(numex,"donnees.dat");
		//printf("%s",numex);
		if (argc == 2)
			optb(argv[1]);
		else
			optb("donnees.dat");
		//printf("\t%d\t%d\t%d\t%d\t",nb_task,nbsepa,nbpas,nbatt1);
		//printf("%ld\t%d\t%d\n",topt,bsup,rsc);
		nnode += nbatt1;
		nsepa += nbsepa;
		nnpas += nbpas;
		temps += topt;
		someh += bsup;
		somop += rsc;
	}
	//printf("\t\t%-7.2lf\t%-7.2lf\t%-7.2lf\t%-7.0lf\t%-7.2lf\t%-7.2lf\n",nsepa/(2-1),
	//nnpas/(2-1),nnode/(2-1),temps/(2-1),someh/(2-1),somop/(2-1));
	printf("%lf\t", temps);
	printf("%d\t", (unsigned int)somop);
	printf("%d\t", LBroot);
	printf("%ld\t", nbsepa);
	//printf("%ld\t", nbPolluted);
	printf("%ld\t", CutActive);
	printf("%ld\n", CutDone);
	//if (CutActive || CutDone)getch();
	printf("cpu:%f\n", toptcpu);
	printf("%lld\n", get_ram_usage());
	if (TimesClean > 0){
		printf("TimesClean:%lld\n", TimesClean);
		printf("NbCleanMinAvgMAx : %lld, %lld, %lld\n", NbCleanMin, NbCleanAvg, NbCleanMax);
	}
	if (SearchStrategy>2)PrintDB();
	printf("NbAddKPerm = %lld\n", NbAddKPerm);
	//! fichier = fopen("PSEDC.txt", "wt"); 
	fichier = fopen("chudc.txt", "wt");
	fprintf(fichier, "%lf\n", toptcpu);
	fprintf(fichier, "%d\n", (unsigned int)somop);
	fprintf(fichier, "%d\n", LBroot);
	fprintf(fichier, "%ld\n", nbsepa);
	fprintf(fichier, "%lf\n", temps);
	fprintf(fichier, "%lld\n", get_ram_usage());
	//fprintf(fichier, "%ld\n",nbPolluted);
	fclose(fichier);

	/*fichier = fopen("statpse1.txt", "wt");
	fprintf(fichier, "%ld\n", CutActive);
	fprintf(fichier, "%ld\n", CutDone);
	fclose(fichier);*/
}


replch(s1, s2)
char *s1, *s2;
{
	char *k, *k1;

	for (k = s2; *k != '\0'; ++k)
		if (isdigit(*k) && !isdigit(*(k - 1))) k1 = k;
	sprintf(s1, "%s", k1);
}


//////////////////////
long long get_ram_usage(){
	static  HANDLE currProc;
	currProc = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS pmc;
	if (K32GetProcessMemoryInfo(currProc, &pmc, sizeof(pmc)))
		return pmc.WorkingSetSize;
	return -1;
}

double get_wall_time(){
	LARGE_INTEGER time, freq;
	if (!QueryPerformanceFrequency(&freq)){
		//  Handle error
		return 0;
	}
	if (!QueryPerformanceCounter(&time)){
		//  Handle error
		return 0;
	}
	return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
	FILETIME a, b, c, d;
	if (GetProcessTimes(GetCurrentProcess(), &a, &b, &c, &d) != 0){
		//  Returns total user time.
		//  Can be tweaked to include kernel times as well.
		double kernelT = (double)(c.dwLowDateTime |
			((unsigned long long)c.dwHighDateTime << 32)) * 0.0000001;
		double userT = (double)(d.dwLowDateTime |
			((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
		//cout << "Kernel = " << kernelT << ";\t User = " << userT<<endl;
		return kernelT + userT;
	}
	else{
		//  Handle error
		fprintf(2, "GetProcessTimes returns 0. Current WallTime = %d\n", get_wall_time());
		return 0;
	}
}