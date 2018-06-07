///////////////////////////////////////////////////////////////////////
//
// Database.h: interface for the Database class. 
// (Working: dont worry of the code in an header file, please!!!!)
//
///////////////////////////////////////////////////////////////////////

#ifndef MYDATABASE

#define MYDATABASE

#include <stdlib.h>
#include <stdio.h>
//#include "beb.h"

#define HASHITEMS 2000	//! initial 500. Can set 2000. The nb of hash positions corresponding to each size of seq. 
						//!HASHITEMS*DBMAXJOBS should be less than MaxDimension.
						//!According to the current hash function, no use to have HASHITEM>sum of id
//#define NB	500		//! only used for declare external var

#define DBMAXJOBS 250	//! initial 500. Can set 200.
#define K_PERM 5 // This parameter is used in DBGenerate() to generate alternative schedules
                 // the K_PERM last jobs are permutated.


#define HASH 1

extern int IsOn; // ON;OFF
struct S_ItemD {

	short Jobs;
	short Seq[DBMAXJOBS];
	int CSum;
	int CMax;
	int Hash;
	char Done;
	char NbUsed;	//! Added for cleaning strategy LUFO
};	

typedef struct S_ItemD ItemD;


extern int Dimension;			// DB Dimension
extern int N_Jobs;				// Number of jobs
extern unsigned nb_task;
extern int 	ed[],pt[];
extern int CutActive, CutDone, CutBef;


extern ItemD *ItemsD;

extern int *Indexes;			 // The indexes vector
	
extern int Starting[HASHITEMS*DBMAXJOBS];	// The starting points for each
											// dimension of the subsequence
											// 0 = empty list.	//! Starting[x] is the hash position for key=x

extern int AllocDB(int Dimmax, int Jobs);
extern int FreeDB(void);

extern signed int DBAdd(int *Seq, int CS, int C2, int Len, int DoNotAdd, int *index2);
							// DoNotAdd =	1 -> only the dominance testing is done.		
							//				-1-> Always add.
							//				0 -> Both.
							// Returns  1 if added (non dominated)
							//			-1 if dominated.
//!MEMO Add (t0, jobset, induced Csum)
signed int DBAddPb(int *Seq, int Len, int t0, int sol, int isLB);

//!MEMO Search a solved pb in Memo. If found, return sol, otherwise return -1.
int DBSearchPb(int *Seq, int Len, int t0, int* isLB);

extern int GetActualDBDimension();
extern int PrintDB(void);

extern int DomTest(int N, int C2a, int CSa, int C2b, int CSb, int minr); // 1  - b is dominated
														// 2  - a is dominated
														// 0  - equal
														// -1 - non dominated


extern int SeqTest(int N, int *A, short *B); // 1 - Same jobs, 0 - no.
	
extern int DBDelete(int *Index, int *LastIndex);	// Delete Items[index], pointed by Lastindex  //! index and lastindex are in valid state at output
											// -1 = first of the list

extern int RemovalCriteria(void);	
						// Removes an item from the db. Used when the DB is full.


extern unsigned int DBGenerate(int *Seq, int CS, int C2, int Len);
// Returns  1 if the input sequence is dominated (the dominating one is added to the DB)
//			0 otherwise.

extern unsigned int DBBuildAndTest(int * Seqtmp, int lenPrefix, int *Set, int lset, int CmaxPrefix, int CSumPermCur);




#endif 
