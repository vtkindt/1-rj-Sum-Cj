#include "database.h"
#include <limits.h>
#include <string.h>
//#include "Stat.h"
//! ===== New =====
#define MEMO //! Enable memo mode. DbAdd will just ... Add
#define NB_MAX	DBMAXJOBS
long long NbAddKPerm = 0;
//! For cleaning
int TimesClean = 0;
extern int CleanStrategy;
extern int DomStrategy;

int IsOn; // ON;OFF
int lastDim = -1;	//!
int CutActive=0, CutDone=0, CutBef=0;

long long NbCleanMin = LLONG_MAX, NbCleanAvg = 0, NbCleanMax = 0;

int Dimension;			// DB Dimension
int N_Jobs;				// Number of jobs

ItemD *ItemsD;

int *Indexes;			 // The indexes vector
	
int Starting[HASHITEMS*DBMAXJOBS];  // The starting points for each //! The starting point for a hash key. A matrix, each line corresponds one N, nb of col is the nb of hashkey positions for entries of size N
									// dimension of the subsequence
									// 0 = beginning of empty list.

int MaxDimension;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int AllocDB(int MaxDim, int Jobs)
{
	int i, RealDimension = MaxDim + 1;  // 1 item is lost for the implementation
										// of the linked list.
    MaxDimension=MaxDim;

	if (MaxDim == 0) {
		IsOn = 0;
		return(0);
	}
	else 
		IsOn = 1;

	N_Jobs = Jobs;
	Dimension = 0;
	printf("Allocating %lf Gb\n", (double)(RealDimension * sizeof(ItemD)) / (double)(1024 * 1024 * 1024));
	ItemsD = (ItemD *)malloc(RealDimension * sizeof(ItemD));

	Indexes = (int *)malloc(RealDimension * sizeof(int));
	
	if (ItemsD == NULL || Indexes == NULL) {
		printf("\nDB_MEMORY_ERROR: Using NO Database\n");
		
		IsOn = 0;
	}

	// Indexes init

	for (i=0;i<MaxDim;i++)
	{
		Indexes[i] = i+1;
		ItemsD[i].NbUsed=ItemsD[i].Done=0;
	}

	Indexes[RealDimension-1] = -1;
	
	for (i=0;i<Jobs*HASHITEMS;i++) /**/
		Starting[i] = -1;

	Starting[0] = 1;

	return(0);
}

int FreeDB()
{
	if (IsOn) {

		free (Indexes);
		free (ItemsD);
	}

	return(0);

}

//////////////////////////////////////////////////////////////////////
// Implementation
//////////////////////////////////////////////////////////////////////
//! 1:A dominates B, 2:B dominates A
int DomTest(int NJ, int C2a, int CSa, int C2b, int CSb, int minr)
{	
	int D, E, n = nb_task - NJ;
	
	// We first test using the condition of Chand et al (1996)

		if (CSa >= CSb) {
			D = n*(C2b - C2a);
			E = CSa - CSb;
			if (D < E) return(2);
		}
		if (CSb >= CSa) {
			D = n*(C2a - C2b);
			E = CSb - CSa;
			if (D < E) return(1);	

		if (C2a == C2b && CSa == CSb) return(0);
		//return(-1);
	}

	if (C2a < C2b && CSa <= CSb) return(1);
	if (C2a <= C2b && CSa < CSb) return(1);
	if (C2a > C2b && CSa >= CSb) return(2);
	if (C2a >= C2b && CSa > CSb) return(2);
	if (C2a == C2b && CSa == CSb) return(0);

    // We next test using the condition of Chu (1992)
    if (CSb<=CSa && n*max(C2a,minr)+CSa>=n*max(C2b,minr)+CSb) return(2);
    if (CSb>=CSa && n*max(C2a,minr)+CSa<=n*max(C2b,minr)+CSb) return(1);


	return(-1);
}

int SeqTest(int NJ, int *A, short *B)
{
	int J[DBMAXJOBS];
	int i, Ret;

	for(i=0;i<N_Jobs;i++) J[i] = 0;

	for(i=0;i<NJ;i++) J[A[i]] = 1;

	Ret = 1;

	for(i=0;i<NJ;i++) 
		if (J[B[i]] == 0) { 
			Ret = 0;
			break;
		}
		
	return(Ret);
}

int DBDelete(int *Index, int *LastIndex)
{
	int k, nj = ItemsD[*Index].Jobs;

	/**/nj = (nj-1) * HASHITEMS + (ItemsD[*Index].Hash % HASHITEMS)+1 ;

	if (*LastIndex == -1) {
	
		k = Starting[0];
		Starting[0] = Starting[nj];
		Starting[nj] = Indexes[Starting[nj]];
		Indexes[Starting[0]] = k;

		*Index = Starting[nj];
		*LastIndex = -1;
	}
	else {

		k = Starting[0];
		Starting[0] = Indexes[*LastIndex];
		Indexes[*LastIndex] = Indexes[*Index];
		Indexes[*Index] = k;

		*Index = Indexes[*LastIndex];
	}

	Dimension--;
	if (Dimension>MaxDimension) printf("\r Database exceeded: %ld",Dimension);
	return(0);
}
int RemovalCriteria()
{
	int i,j;
	long long totalNbDeleted = 0;
	if (CleanStrategy == 0){//! default cleaning implemented on 2004, the head of the list of the fewest jobs
		for (i = 1; i < N_Jobs*HASHITEMS; i++) /**/
			if (Starting[i] != -1) break;

		// Deleting the first item of the list with less jobs.

		i = Starting[i];
		j = -1;
		DBDelete(&i, &j);
		totalNbDeleted++;
	}
	else if (CleanStrategy == 1 || CleanStrategy==3)//! Least used first
	{
		//! Remove all unused entries and decrease the ctr by one
		for (i = 1; i < N_Jobs*HASHITEMS; i++) {/**/
			if (Starting[i] != -1)//! a list of items found
			{
				int curr = Starting[i], last = -1;
				while (curr != -1){
					if (ItemsD[curr].NbUsed <= 0)
					{
						DBDelete(&curr, &last);
						totalNbDeleted++;
					}
					else{
						ItemsD[curr].NbUsed--;
						last = curr;
						curr = Indexes[curr];
					}
				}
			}
		}
		if (totalNbDeleted <= 0){
			// Deleting the last list with most jobs.
			for (i = N_Jobs*HASHITEMS - 1; i >= 1; i--) /**/
				if (Starting[i] != -1) break;
			int curr = Starting[i];
			int last = -1;
			while (curr != -1){
				DBDelete(&curr, &last);
				totalNbDeleted++;
			}
		}
	}
	else if (CleanStrategy == 2 )//! Longest first
	{
		//! Search from end to begin x positions until find a non empty list, then clean (at least) x items from then. 
		int x = 0;
		for (i = N_Jobs*HASHITEMS - 1; i>=1 ; i--) {/**/
			if (Starting[i] != -1)//! a list of items found
			{
				if(x==0)x = N_Jobs*HASHITEMS - i;//! nb of steps performed to find a list.
				int curr = Starting[i], last = -1;
				while (curr != -1){
					DBDelete(&curr, &last);
					totalNbDeleted++;
				}
				if (totalNbDeleted > x)break;
			}
		}
	}
	TimesClean++;
	if (TimesClean % 1 == 0)printf("TimesClean=%d, NbCleanMin=%d, NbCleanAvg=%d, NbCleanMax=%d\n", TimesClean, NbCleanMin, NbCleanAvg, NbCleanMax);
	if (NbCleanMax < totalNbDeleted)NbCleanMax = totalNbDeleted;
	if (NbCleanMin > totalNbDeleted)NbCleanMin = totalNbDeleted;
	NbCleanAvg = NbCleanAvg + (totalNbDeleted - NbCleanAvg)*1.0 / TimesClean;
	return(Starting[0]); // Returns the freed item.
}


void DBSearch(int *Seq, int CS, int C2, int Len, int *Rindex, int *lastindex)
{
	int i, index, l_index, Sum;
	
	int NH_NJobs;
	int N_NJobs = Len;
	int N_CMax  = C2; 
	int N_CSum  = CS;

	if (IsOn == 0) return;	

	// Search for a dominated sequence to be deleted
	// or for a dominating sequence.
	
	Sum = 0;	
	if (HASH) {
		for (i=0; i<N_NJobs; i++) Sum += Seq[i];
	}

	NH_NJobs = (N_NJobs-1)*HASHITEMS + (Sum % HASHITEMS)+1; //! starting index for seq of size Njobs + offset given by key%HASHITEMS

	index = Starting[NH_NJobs]; /**/
	l_index = -1;

	while(index != -1) {
	
			if (Sum == ItemsD[index].Hash) {

				if ((N_CMax==ItemsD[index].CMax)&&(N_CSum==ItemsD[index].CSum)
					&&(SeqTest(N_NJobs, Seq, ItemsD[index].Seq)))
				{
				 *Rindex=index;
				 *lastindex=l_index;
				 return;
				}
			}
		
			l_index = index;
			index = Indexes[index];
		}
 *Rindex=index;
 *lastindex=l_index;
}

//! Search a solved pb in Memo. If found, return sol, otherwise return -1
int DBSearchPb(int *Seq, int Len, int t0, int * isLB)
{
	int i, index, l_index, Sum, nbUsedNewItem = 0;

	int NH_NJobs;
	int N_NJobs = Len;

	if (IsOn == 0) return -1;

	Sum = 0;
	if (HASH) {
		for (i = 0; i < N_NJobs; i++) Sum += Seq[i];
	}

	NH_NJobs = (N_NJobs - 1)*HASHITEMS + (Sum % HASHITEMS) + 1; //! starting index for seq of size Njobs + offset given by key%HASHITEMS

	index = Starting[NH_NJobs]; /**/
	l_index = -1;

	while (index != -1) {

		if (Sum == ItemsD[index].Hash) {
			//! Cmax is used to store t0
			if ((t0 == ItemsD[index].CMax) && (SetTest(N_NJobs, Seq, ItemsD[index].Seq)==1))
			{
				*isLB = ItemsD[index].Done;
				++ ItemsD[index].NbUsed;
				return  ItemsD[index].CSum;
			}
		}
		index = Indexes[index];
	}
	return -1;
}

signed int DBAdd(int *Seq, int CS, int C2, int Len, int DoNotAdd, int *index2)
{
	int i, index, l_index, Sum, nbUsedNewItem=0;
	int NH_NJobs;
	int N_NJobs = Len;
	int N_CMax  = C2; 
	int N_CSum  = CS;
	int minr;

	if (IsOn == 0) return(1);	

	// Search for a dominated sequence to be deleted
	// or for a dominating sequence.
	
	Sum = 0;	
	if (HASH) {
		for (i=0; i<N_NJobs; i++) Sum += Seq[i];
	}

	NH_NJobs = (N_NJobs-1)*HASHITEMS + (Sum % HASHITEMS)+1;

	if (DoNotAdd != -1) {

		index = Starting[NH_NJobs]; /**/
		l_index = -1;
		minr=999999;
		for (i=Len;i<nb_task;i++) if (ed[Seq[i]]<minr) minr=ed[Seq[i]];

		while(index != -1) {
	
			if (Sum == ItemsD[index].Hash) {
				i = DomTest(N_NJobs, N_CMax, N_CSum, ItemsD[index].CMax, ItemsD[index].CSum,minr);

				if ((i == 2)&&((*index2)==0)&&(ItemsD[index].Done!=100)) {
					if ((SeqTest(N_NJobs, Seq, ItemsD[index].Seq)) &&(ItemsD[index].Jobs==N_NJobs))	/*Cut by active*/
					{
						*index2=index;
						++ ItemsD[index].NbUsed;//! Added for LUFO
						return(-1); // N is Dominated by a seq. in the DB
					}
				}
				if ((i == 2)&&((*index2)==100)&&(ItemsD[index].Done==100)) {
					if ((SeqTest(N_NJobs, Seq, ItemsD[index].Seq)) &&(ItemsD[index].Jobs==N_NJobs))	/*Cut by done*/
					{
						*index2=index;
						++ ItemsD[index].NbUsed;
						return(-1); // N is Dominated by a seq. in the DB
					}
				}

				if (i == 1) {
					if ((SeqTest(N_NJobs, Seq, ItemsD[index].Seq)) &&(ItemsD[index].Jobs==N_NJobs)) { /**/
						if (CleanStrategy == 3)nbUsedNewItem += ItemsD[index].NbUsed;
						DBDelete(&index, &l_index); // Delete index from the db
						continue;					// and move the indexes.
					}
				}
				if (DomStrategy & 4){
					if (i == 0 && N_CMax == ItemsD[index].CMax && N_CSum == ItemsD[index].CSum){
						// This test is added by VTkindt (04/03/2017) to cope with the case where the current solutioh has already been added to the database by DBGenerate()
						if ((SeqTest(N_NJobs, Seq, ItemsD[index].Seq)) && (ItemsD[index].Jobs == N_NJobs)) { /**/

							DBDelete(&index, &l_index); // Delete index from the db
							continue;					// and move the indexes.
						}
					}
				}
			}
		
			l_index = index;
			index = Indexes[index];
		}

	}

	// Free list insertion (if possible)

	index = Starting[0];	//! The position where to insert

	if (index == -1) {
		index = RemovalCriteria(); // db full
		static char msg = 0;
		if (msg==0) {
			msg = 1;
			printf("DB full. RAM=%lld", get_ram_usage());
		}
	}
	
	if (DoNotAdd == 1) return(1);

	// REAL INSERTION
	
	// New empty list starting
	Starting[0] = Indexes[Starting[0]];
		
	// New item.jobs startings

	i = Starting[NH_NJobs];			//! i=beginning of list
	Starting[NH_NJobs] = index;		//! beginning list = new position
	Indexes[index] = i;				//! new position -> old beginning of list => inserted at the beginnging
	
	// Copying data into Items[index]
		
	ItemsD[index].Jobs = N_NJobs;

	for (i=0; i<N_NJobs; i++) ItemsD[index].Seq[i] = Seq[i]; /**/

	ItemsD[index].Hash = Sum;
	ItemsD[index].CSum = N_CSum;
	ItemsD[index].CMax = N_CMax;
	ItemsD[index].Done = 100;
	ItemsD[index].NbUsed = nbUsedNewItem;

	Dimension++;
	//!
	if (lastDim<Dimension && Dimension % 100000 == 0){
		lastDim = Dimension;
		printf("Current DBDimension : %d\n", Dimension);
	}

	return(1);
}

//! Return 1 if A and B are the same set
int SetTest(int NJ, int *A, short *B)
{
	int J[DBMAXJOBS+1];
	int i, Ret;

	for (i = 0; i < N_Jobs+1; i++) J[i] = 0;

	for (i = 0; i < NJ; i++) J[A[i]] = 1;

	Ret = 1;

	for (i = 0; i < NJ; i++)
		if (J[B[i]] == 0) {
			Ret = 0;
			break;
		}
	return(Ret);
}

//! Add (t0, jobset, induced Csum)
signed int DBAddPb(int *Seq, int Len, int t0, int sol, int isLB){
	int i, index, l_index, Sum;


	int NH_NJobs;
	int N_NJobs = Len;

	if (IsOn == 0) return(1);
	Sum = 0;
	if (HASH) {
		for (i = 0; i < N_NJobs; i++) Sum += Seq[i];
	}

	// Free list insertion (if possible)
	index = Starting[0];	//! The position where to insert
	if (index == -1)
	{
		//printf("Database is full\n");
		index = RemovalCriteria(); // db full
	}
	// New empty list starting
	Starting[0] = Indexes[Starting[0]];
	// REAL INSERTION
	NH_NJobs = (N_NJobs - 1)*HASHITEMS + (Sum % HASHITEMS) + 1;
	// New item.jobs startings
	i = Starting[NH_NJobs];			//! i=beginning of list
	Starting[NH_NJobs] = index;		//! beginning list = new position
	Indexes[index] = i;				//! new position -> old beginning of list => inserted at the beginnging

	// Copying data into Items[index]
	ItemsD[index].Jobs = N_NJobs;

	for (i = 0; i < N_NJobs; i++) ItemsD[index].Seq[i] = Seq[i]; /**/

	ItemsD[index].Hash = Sum;
	ItemsD[index].CMax = t0;	//! Use Cmax as t0
	ItemsD[index].CSum = sol;	//! Use Csum as sol value
	ItemsD[index].Done = isLB;  //! Use Done to indicate in sol memo whether the value is a lb (or the optimal)
	ItemsD[index].NbUsed = 0;

	Dimension++;
	//!
	if (lastDim < Dimension && Dimension % 100000 == 0){
		lastDim = Dimension;
		printf("Current DBDimension : %d\n", Dimension);
	}

	return(1);
}

int GetActualDBDimension()
{
	return(Dimension);
}

int PrintDB(void)
{
	int i,j;
	
	printf("Empty Slots (jamais servis) = ");
	
	i = Starting[0];
	j = 0;

	while(i != -1) {
		j++;
		i = Indexes[i];
	}
	
	printf("%d\nDimension=%d\n", j, Dimension);
	return 0;
	for(i=1; i<N_Jobs; i++)//! Pb: the indices of Starting are hashed values (max:N_job*HASHITEMS) instead of jobsize.
		
		if (Starting[i] != -1) {
		
			printf("%d Jobs:\n", i);

			j = Starting[i];

			while(j != -1) {

				printf("  C2 = %3d, CS = %3d \n", ItemsD[j].CMax, ItemsD[j].CSum);
			
			
				j = Indexes[j];
			}
		}

	printf("=================================================\n");
	return(0);
}

//! TODO LEi: we should somehow avoid generating on the same jobset
int seqPermBest[K_PERM], CmaxPrefix0 = 0, CSum0 = 0, Cmax0 = 0, CSumPerm0 = 0;//! currently best records,  //! seq: the best seq among all k perm
int lenPrefix0 = 0, lenPerm0 = 0, len0 = 0, minr0;							  //! global invariants
unsigned int DBGenerate(int *Seq, int CS, int C2, int Len)
{
	int i;
	//int *Seqtmp, *Set;
	int Seqtmp[NB_MAX], Set[K_PERM];	//! tmp seq during generation & the last k jobs

	CSum0 = CS;
	CmaxPrefix0 = 0, lenPrefix0 = 0, lenPerm0=0;
	len0 = Len;
	minr0 = 999999;

	if (IsOn == 0) return(1);
	// Search for a sequence dominating the input sequence
	for (i = Len; i<nb_task; i++) if (ed[Seq[i]]<minr0) minr0 = ed[Seq[i]];

	// Phase 1: studying permutations of the k last jobs in Seq
	//Seqtmp = (int *)malloc(sizeof(int)*Len);
	//Set = (int *)malloc(sizeof(int)*K_PERM);
	lenPrefix0 = (signed int)Len - (signed int)K_PERM;
	if (lenPrefix0 < 0)lenPrefix0 = 0;
	for (i = 0; i < lenPrefix0; i++){
		Seqtmp[i] = Seq[i];
		CmaxPrefix0 = max(CmaxPrefix0, ed[Seq[i]]) + pt[Seq[i]];
	}
	lenPerm0 = Len - lenPrefix0;
	for (i = 0; i < lenPerm0; i++){
		seqPermBest[i] = Set[i] = Seq[lenPrefix0 + i];
	}
	CSumPerm0 = 0;
	Cmax0 = CmaxPrefix0;
	for (i=lenPrefix0; i < Len; ++i){
		Cmax0 = max(Cmax0, ed[Seq[i]]) + pt[Seq[i]];
		CSumPerm0 += Cmax0;
	}
	if (Cmax0 != C2){ printf("Cmax0 != N_CMax!\n"); exit(1); }

	unsigned int Status = DBBuildAndTest(Seqtmp, lenPrefix0, Set, lenPerm0, CmaxPrefix0, 0);
	if (Status){
		NbAddKPerm++;
		if ((DomStrategy & 10) == 10){// add the best to db.
			for (i = 0; i < lenPerm0; ++i){
				Seqtmp[lenPrefix0+i] = seqPermBest[i];
				DBAdd(Seqtmp, CSum0 , Cmax0, len0, -1, NULL);
			}
		}
	}
	//free(Seqtmp);
	//free(Set);
	return (Status);
}


//! Len:lSR-lSet
unsigned int DBBuildAndTest(int * Seqtmp, int lenPrefix, int *Set, int lset, int CmaxPrefix, int CSumPermCur)
{
	unsigned int uiLoop, uiLoop2, Status=0;
	int Settmp[K_PERM];
	int lsettmp;
	int Cmax, i;

	if (lset > 0)
	{
		//Settmp = (int *)malloc(sizeof(int)*(lset - 1));
		for (uiLoop = 0; uiLoop < lset; uiLoop++)
		{
			Seqtmp[lenPrefix] = Set[uiLoop];//! a new job fixed at the first position
			lsettmp = 0;
			for (uiLoop2 = 0; uiLoop2 < lset; uiLoop2++)
				if (uiLoop != uiLoop2)
				{
					Settmp[lsettmp] = Set[uiLoop2];
					lsettmp++;
				}
			int CmaxPrefixNew = max(CmaxPrefix, ed[Set[uiLoop]]) + pt[Set[uiLoop]];
			Status = DBBuildAndTest(Seqtmp, lenPrefix + 1, Settmp, lsettmp, CmaxPrefixNew, CSumPermCur + CmaxPrefixNew)
				|| Status;
			if ( !(DomStrategy & 2) && Status == 1)
			{// If at least one improving sequence is found we stop exploring the neighborhood
				//free(Settmp);
				return Status; 
			}
		
		}
		//free(Settmp);
		if(DomStrategy & 2) 
			return Status;
		else return 0;
	}
	else
	{ 
		// We need to evaluate if the sequence Seqtmp improves upon seqBestPerm on the last kperm positions. The reference: seqBestPerm, Cmax0, CsumPerm0
		if (lenPrefix != len0)printf("Problem!\n");
		Cmax = CmaxPrefix;
		//CsumPerm = 0;
		//for (uiLoop = lenPrefix; uiLoop < len0; uiLoop++)
		//{
		//	if (CmaxPerm < ed[Seqtmp[uiLoop]]) CmaxPerm = ed[Seqtmp[uiLoop]] + pt[Seqtmp[uiLoop]];
		//	else CmaxPerm += pt[Seqtmp[uiLoop]];
		//	CsumPerm += CmaxPerm;
		//}

		Status = DomTest(len0, Cmax0, CSum0, Cmax, CSum0 - CSumPerm0 + CSumPermCur, minr0);
		if (Status == 2)
		{   // The reference sequence is dominated by the built sequence Seqtmp
			// if DomStra&0x10 == 0, the latter is added to the DB and we are done
			if (!(DomStrategy & 2) && (DomStrategy & 8)){
				DBAdd(Seqtmp, CSum0 - CSumPerm0 + CSumPermCur, Cmax, len0, -1, NULL);
			}
			// else, we update the current best perm and we add the best when all perm are enumerated
			else{
				for (int i = 0; i < lenPerm0; ++i)
					seqPermBest[i] = Seqtmp[lenPrefix0 + i];
				Cmax0 = Cmax;
				CSum0 = CSum0 - CSumPerm0 + CSumPermCur;
				CSumPerm0 = CSumPermCur;
			}
			return(1);
		}
		else
		{
			return(0);
		}
	}
}