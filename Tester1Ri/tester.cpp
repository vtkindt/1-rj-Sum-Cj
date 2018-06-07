#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include "testutils.h"
#pragma warning(disable:4996)
#define PMAX 100

// Important : init par.rameters and stats
Config c;
Stat s;

// This function randomly generates one instance and copy it to /data_gen
long int generer(int nbtrav, double range, char* copyto);

// arg1=strategy; arg2=startsize; arg3=exe path
void main(int argn, char*argv[])
{
	//argn = 2; char *argvv[] = { "x", "tester_config.txt" }; argv = argvv;
	if (argn == 2) {
		c.readConfig(argv[1]);
		c.printConfig();
	}else if (argn != 6){ 
		printf("Usage by config file : %s [configfile]\n", *argv);
		printf("Usage by cmd args    : %s [BrStra] [StartSize] [exe] [CleaningStra] [DomStra]\n", *argv);
		printf("Usage for sipsi      : %s 10 [StartSize] [exe] 0 0\n", *argv);
		exit(0);
	}
	else {
		c.VERSION = atoi(argv[1]);
		c.SIZE_START = atoi(argv[2]);
		c.SOLVER_PATH=string(argv[3]);
		c.CLEAN_STRA = atoi(argv[4]);
		c.DOM_STRA = atoi(argv[5]);
	}

	setvbuf(stdout, NULL, _IONBF, 0);// No buffering
	
	// Print headers to files
	FILE *fichier = fopen(c.OUTPUT_FILE.c_str(), "at");
	fprintf(fichier, "Tests of Sipsi\n Done on a PC HPZ400 3.06Ghz 8G RAM\n");
	// Header info should always be updated.
	if(c.VERSION==10)
		fprintf(fichier, "Ins\t Sol\t TCpu\t #Nodes\t #Arcs\t TWall\t Ram\n");
	else 
		fprintf(fichier, "Ins\t Sol\t TCpu\t #Nodes\t TWall\t Ram\n");
	fprintf(fichier, "--------------------------------------------------------------------\n");
	fclose(fichier);

	fichier = fopen(c.STAT_FILE.c_str(), "at");
	fprintf(fichier, "These tests are done on a PC HPZ400 3.06Ghz 8G RAM\n");
	fprintf(fichier, "Number of iterations for a fixed value of n : %d\n\n", c.N_INS_PER_SETTING);
	fprintf(fichier, "Ins  nbIns  Nodes_min   Nodes_moy  Nodes_max  Time_min  Time_moy  Time_max\n");
	fprintf(fichier, "--------------------------------------------------------------------\n");
	fclose(fichier);

	FILE *fResIn;
	if (c.VERSION == 30) {
		fResIn = fopen(c.INPUT_FILE.c_str(), "r");
		if (!fResIn) {
			fprintf(stderr, "Problem on opening input file %s. Exiting...\n", c.INPUT_FILE.c_str());
			exit(0);
		}
	}
	char datadir[50], buf[50];
	double ranges[] = { 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 3.0 };
	const char* strRanges[] = { "0.2", "0.4", "0.6", "0.8", "1.0", "1.25", "1.5", "1.75", "2.0", "3.0" };
	for (int insSize = c.SIZE_START; insSize <= c.MAX_SIZE_PB; insSize += c.SIZE_STEP)
	{
		int insId = 0;
		if (c.VERSION == 20) {
			sprintf(datadir, "md data_gen\\%d\\", insSize);
			system(datadir);
		}
		//srand(1);
		for (int indRange = 0; indRange < 10; ++indRange)
		{
			double range = ranges[indRange];
			printf("N,R : %ld,%lf\n", insSize, range);
			for (int j = 0; j<c.N_INS_PER_SETTING; j++)
			{
				//printf("Jeux n°%d\n", j + 1);
				++insId;
				if (insSize == c.SIZE_START && insId < c.INS_START)continue;
				//printf("data\\%d\\%d.txt", insSize, insId);
				int indInRange = (insId -1) % c.N_INS_PER_SETTING+1;
				if (c.VERSION == 20) {
					sprintf_s(datadir, "data_gen\\%d\\SDT_%d_%s_%d.txt", insSize, insSize,strRanges[indRange], indInRange);
				}else
					sprintf_s(datadir, "data\\%d\\SDT_%d_%s_%d.txt", insSize, insSize, strRanges[indRange], indInRange);

				generer(insSize, range, datadir);
				if (c.VERSION == 20)continue;	// Data generation

				Result r;						
				if (c.VERSION == 30) {
					//read from res file
					int szTmp=insSize-1, idTmp=insId-1;
					int ctr = 0;
					while (true) {
						//if(c.NB_COL_INPUT_FILE == 8)
							ctr=fscanf(fResIn, "data\\%d\\%d.txt\t%d\t%lf\t%lld\t%lld\t%lf\t%lld\n",
								&szTmp, &idTmp, &r.sol, &r.tCpu, &r.nbNodes, &r.nbArcs, &r.tWall, &r.ram);
						printf("%d, %d, %d\n", ctr, szTmp, idTmp);
						if (ctr!=8 || (szTmp == insSize && idTmp == insId)) break;
						else if (szTmp < insSize || (szTmp==insSize && idTmp<insId)) continue;
						else {
							fprintf(stderr,"Cannot find result from file for instance %s. Exiting...\n", datadir);
							exit(0);
						}
					}
					if (ctr != 8) { indRange = 10; c.MAX_SIZE_PB = insSize - 1; --insId; break; }
					s.addResult(r, insSize);
				}
				else if (c.VERSION == 10) { // Sipsi
					sprintf_s(buf, "%d", insSize);
					_spawnl(P_WAIT, c.SOLVER_PATH.c_str(), c.SOLVER_PATH.c_str(), datadir, buf,NULL);
					fichier = fopen(c.SOLVER_SOL_FILE.c_str(), "rt");
					fscanf_s(fichier, "%lf\n",  &r.tCpu);
					fscanf_s(fichier, "%d\n", &r.sol);
					//fscanf_s(fichier, "%d\n", &lb);
					fscanf_s(fichier, "%lld\n", &r.nbNodes);
					fscanf_s(fichier, "%lld\n", &r.nbArcs);
					fscanf_s(fichier, "%lf\n",  &r.tWall);
					fscanf_s(fichier, "%lld\n", &r.ram);
					fscanf_s(fichier, "%d\n",&r.retCode);
					fclose(fichier);
					s.addResult(r, insSize);

					fichier = fopen(c.OUTPUT_FILE.c_str(), "at");
					fprintf(fichier, "%s\t%d\t%lf\t%lld\t%lld\t%lf\t%lld\n", 
						datadir, r.sol, r.tCpu, r.nbNodes, r.nbArcs, r.tWall, r.ram);
					fclose(fichier);
					//! Stop progr.ram
					if (r.tCpu > c.TIME_LIM) { indRange = 10; c.MAX_SIZE_PB = insSize - 1; break; }//exiting
				}
				// Exécution de la PSE de CHU
				else //if (c.VERSION < 10)
				{
					//printf("Chudc is running...\n");
					fichier = fopen("chudc.ini", "wt");
					fprintf(fichier, "%d\n", c.VERSION);
					fprintf(fichier, "%d\n", c.CLEAN_STRA);
					fprintf(fichier, "%d\n", c.DOM_STRA);
					fclose(fichier);
					_spawnl(P_WAIT, c.SOLVER_PATH.c_str(), c.SOLVER_PATH.c_str(), datadir, NULL);
					fichier = fopen(c.SOLVER_SOL_FILE.c_str(), "rt");
					fscanf_s(fichier, "%lf\n", &r.tCpu);
					fscanf_s(fichier, "%d\n", &r.sol);
					fscanf_s(fichier, "%d\n", &r.lb);
					fscanf_s(fichier, "%lld\n", &r.nbNodes);
					fscanf_s(fichier, "%lf\n", &r.tWall);
					fscanf_s(fichier, "%lld\n", &r.ram);
					//fscanf_s(fichier, "%ld\n", &r.nbPolluted);
					fclose(fichier);
					s.addResult(r, insSize);

					fichier = fopen(c.OUTPUT_FILE.c_str(), "at");
					fprintf(fichier, "%s\t%d\t%.2f\t%lld\t%.2f\t%lld\n", 
						datadir, r.sol, r.tCpu, r.nbNodes, r.tWall, r.ram);
					fclose(fichier);
					//! Stop progr.ram
					if (r.tCpu > c.TIME_LIM){ indRange = 10; c.MAX_SIZE_PB = insSize - 1; break; }
				}
			}//for iteration
		}//one size finished

		//! Results by size
		auto minAvgMaxNodes = s.getMinAvgMax(insSize, &Result::nbNodes);
		auto minAvgMaxT = s.getMinAvgMax(insSize, &Result::tCpu);
		fichier = fopen(c.STAT_FILE.c_str(), "at");
		fprintf(fichier, "%d\t%d\t%lld\t%.2lf\t%lld\t%.3lf\t%.3lf\t%.3lf\n", 
			insSize, insId, get<0>(minAvgMaxNodes), get<1>(minAvgMaxNodes), get<2>(minAvgMaxNodes), 
			get<0>(minAvgMaxT), get<1>(minAvgMaxT), get<2>(minAvgMaxT));
		fclose(fichier);
	}
}

// ========
// This function randomly generates one instance and copy it to /data_gen
long int generer(int nbtrav, double range, char* copyto)
{
	if (c.VERSION != 20) return -1;

	int i;
	int ri, som = 0;
	int * pi = new int[c.MAX_SIZE_PB];
	FILE *fichier, *fcopy;

	fichier = fopen(c.INS_FILE.c_str(), "wt");
	if (copyto != NULL)
	{
		fcopy = fopen(copyto, "wt");
	}
	int sumri = 0;
	for (i = 0; i<nbtrav; i++)
	{
		pi[i] = int(((float)rand() / (float)RAND_MAX) * (PMAX-1) + 1);
		som += pi[i];
	}
	for (i = 0; i<nbtrav; i++)
	{
		ri = int(((float)rand() / (float)RAND_MAX)*(PMAX+1)*0.5*nbtrav*range);//! the old 50.5 is replaced by (pmax+1)/2, by *n, representing the expected total processing time
		sumri += ri;
		fprintf(fichier, "%ld %ld %ld\n", i + 1, ri, pi[i]);
		if (copyto != NULL)
			fprintf(fcopy, "%ld %ld %ld\n", i + 1, ri, pi[i]);
	}
	fclose(fichier);
	if (copyto != NULL) fclose(fcopy);
	printf("sumri=%ld\n", sumri);
	delete[] pi;
	return(sumri);
}