# Organisation

**AllProjects**: open the vs solution containing all 4 projects.

**Exact**: code project of the BB

**LBsrptImp**: code project of the *srptImp* algorithm

**sipsi**: code project of Tanaka's algorithm

**Tester1Ri**: code project of the tester. (the solution file is inside)

**Output**: the exe file of all projects are copied to this folder. Use this folder to run any program. The dataset is also there. **Extract the data first**.

# Parameters

### To run the tester:

Put parameters in **tester_config.txt**, then run **tester.exe  tester_config.txt**

The meaning of parameters can be found in **tester_config.txt**. 

Or in **Tester1Ri/testutils.h**. Check the **Config** class.

The **data.zip**  in **Output** folder should be decompressed first.

If need to regenerate data, verify the macro **PMAX** in *tester.cpp*.

### To run the BB (chudc.exe) separately:

- "chudc.exe [file=donnees.dat]"

- the chudc.ini must exist. One parameter per line

- Line 1: Search Strategy

   0, classical depth first
   1, classical best first
   2, classical breadth first
   3, depth first memo
   4, best first memo
   5, breadth first memo
   6, solution/lb memo (depth first)

- Line 2: Clean Strategy. Put 3 to use LUFO.

- Line 3: Dominance Strategy. Put 0 to disable or 12 to enable predictive dominance (by checking k-permutation)

### To run the sipsi algorithm (sipsi.exe) separately:

- "sipsi.exe input_file_name nbJobs"

  