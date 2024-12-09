Overview - K-means Clustering Implementation

Project Authors:
*Sophie Prahia
*Romy Shufman

This project was created as part of a university software project course. It implements the K-means clustering algorithm in Python. Developed in 2022, it provides a robust solution for unsupervised machine learning and data segmentation.
Unsupervised machine learning clustering algorithm, Supports multiple distance metrics, Handles multi-dimensional data, Customizable number of clusters

├── spkmeans.py      # Python interface
├── spkmeans.c       # C implementation
├── spkmeans.h       # C header file
├── spkmeansmodule.c # Python C API wrapper
├── setup.py         # Build configuration
└── Makefile        # Compilation script

spkmeans.py (Python Interface) - This is the main Python script that users interact with. It handles input/output, data processing, and calling the C functions.It Uses import mykmeanssp to access the C functions.

spkmeans.c and spkmeans.h (C Implementation)- 
spkmeans.c: Contains the actual implementation of the algorithms (wam, ddg, jacobi, etc.)
spkmeans.h: Declares the functions and structures that spkmeans.c implements


spkmeansmodule.c (The "Bridge")- This is the *translator* between Python and C. It contains wrapper functions that:
1. Take Python data → Convert to C format
2. Call the C functions
3.Take C results → Convert back to Python format


setup.py (Building Instructions)- This file tells Python how to build your C code into something it can import. Specifies which files to compile and link together. When you run python setup.py build, it creates the importable module

Makefile (Compilation Helper)- Helps compile and manage the C code components. Makes it easier to build and clean the project

The Flow:
CopyUser → spkmeans.py → mykmeanssp → spkmeansmodule.c → spkmeans.c → Results
      (Python)       (Import)      ("Bridge")        (C code)

Python loads spkmeans.py
When it sees import mykmeanssp, it loads the *compiled* C module
When you call a function like mykmeanssp.wam():

The call goes through spkmeansmodule.c,
Gets translated to C format,
Runs your C implementation,
Results come back through the translation layer,
Python gets the results

Algorithm Details:

Input: Dataset, number of clusters (K)
Initialize K centroids
Assign points to nearest centroid
Recalculate centroids
Repeat until convergence





gcc spkmeans.c -o spkmeans -lm
spkmeans wam text.txt
spkmeans ddg text.txt
spkmeans jacobi jacobi1.txt
python3 setup.py build_ext --inplace
python3 spkmeans.py 3 spk text.txt

python3 spkmeans.py wam testbatch/test_batch/test1.txt > testbatch/ans_batch/test1_wam.txt

python3 spkmeans.py ddg testbatch/test_batch/test10.txt > testbatch/ans_batch/test10_ddg.txt

python3 spkmeans.py gl testbatch/test_batch/test10.txt > testbatch/ans_batch/test10_gl.txt

python3 spkmeans.py jacobi testbatch/test_batch/test1_j.txt > testbatch/ans_batch/test1_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test2_j.txt > testbatch/ans_batch/test2_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test3_j.txt > testbatch/ans_batch/test3_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test4_j.txt > testbatch/ans_batch/test4_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test5_j.txt > testbatch/ans_batch/test5_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test6_j.txt > testbatch/ans_batch/test6_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test7_j.txt > testbatch/ans_batch/test7_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test8_j.txt > testbatch/ans_batch/test8_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test9_j.txt > testbatch/ans_batch/test9_j_ans.txt
python3 spkmeans.py jacobi testbatch/test_batch/test10_j.txt > testbatch/ans_batch/test10_j_ans.txt
