Overview - K-means Clustering Implementation

This project implements the K-means clustering algorithm in Python. Developed in 2021, it provides a robust solution for unsupervised machine learning and data segmentation.
Features

Unsupervised machine learning clustering algorithm
Supports multiple distance metrics
Handles multi-dimensional data
Visualizes clustering results
Customizable number of clusters

Key Components:

K-means clustering algorithm
Distance calculation methods
Centroid initialization strategies
Convergence detection
Visualization utilities

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
