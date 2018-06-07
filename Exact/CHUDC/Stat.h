#ifndef STAT_H
#define STAT_H
#include <iostream>
#include <map>
#include <vector>
#include <Windows.h>
#include "psapi.h"
using namespace std;

template <typename T = long long>
class Stat{   // Stat for one criterion
public:
	T mi, ma, sum;
	double av;
	long long maxRam;
	double lastCpuTime, lastWallTime;
public:
	Stat(){ lastWallTime = lastCpuTime = maxRam = -1; }
	vector<T> data;
	int add(T v){
		data.push_back(v);
		if (data.size() == 1){ mi = av = ma = sum = data.front(); return 0; }
		sum += data.back();
		av = 1.0 * sum / data.size();
		if (data.back() < mi){ mi = data.back(); return -1; }
		if (data.back() > ma){ ma = data.back(); return 1; }
	}
	// Get min, avg, max
	T getMin(){ if (data.empty()) throw out_of_range("No data no min."); else return mi; }
	double getAvg(){ if (data.empty()) throw out_of_range("No data no avg."); else return av; }
	T getMax(){ if (data.empty()) throw out_of_range("No data no max."); else return ma; }

	// Methods for querying time/ram
	static long long getRam(){
		static  HANDLE currProc;
		currProc = GetCurrentProcess();
		PROCESS_MEMORY_COUNTERS pmc;
		if (K32GetProcessMemoryInfo(currProc, &pmc, sizeof(pmc)))
			return pmc.WorkingSetSize;
		throw runtime_error("getRam->K32GetProcessMemoryInfo return false.");
		return -1;
	}
	static double getWallTime(){
		LARGE_INTEGER time, freq;
		if (!QueryPerformanceFrequency(&freq)){
			//  Handle error
			throw runtime_error("getWallTime->QueryPerformanceFrequency return false.");
			return 0;
		}
		if (!QueryPerformanceCounter(&time)){
			//  Handle error
			throw runtime_error("getWllTime->QueryPerformanceCounter return false.");
			return 0;
		}
		return (double)time.QuadPart / freq.QuadPart;
	}
	static double getCpuTime(){
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
			fprintf(stderr, "GetProcessTimes returns 0. Current WallTime = %d\n", getWallTime());
			return 0;
		}
	}
	double getWallDuration(){
		if (lastWallTime == -1){ lastWallTime = getWallTime(); return 0; }
		else { double res = lastWallTime; lastWallTime = getWallTime(); return lastWallTime - res; }
	}
	double getCpuDuration(){
		if (lastCpuTime == -1){ lastCpuTime = getCpuTime(); return 0; }
		else { double res = lastCpuTime; lastCpuTime = getCpuTime(); return lastCpuTime - res; }
	}
	long long updateMaxRam(){
		long long ram = getRam();
		if (ram > maxRam) maxRam = ram;
		return ram;
	}
};


#endif