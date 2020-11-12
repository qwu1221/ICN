# ICN
This is an implementation of ICN detection for gene co-expression network (GCN) in Matlab. 
In constrast to the community structure, we consider a pair of communities be interconnected if a subset of genes from one community is correlated with a subset of genes from another community. Please refer to the paper "Extracting Interconnected Communities in Gene co-expression Networks" for details. 

# Codes
The algorithm includes 3 steps. 

In step 1, we detect a set of densely connected, non-overlapping communities as the backbone of ICN structure. The implementation refers to the "NICE.m" in a previous package https://github.com/qwu1221/ICN/blob/master/README.md.

In step 2, we test whether a pair of communities are interconnected based on KL divergence. The testing procedure is described in algoirthm 1 of the paper, and the implementation is in "KLtest.m".

In step 3: we identify connecting edges for a pair of interconnected communities based on the results from step 2. We implement the algorithm using "InterCut.m";

# Examples
We include a small simulation example with ICN detection "demo_3block.m". This file includes the generating mechanism of a simulation data set and the implementation of all 3 steps. To see the demonstrative document with codes and output, please refer to "demo_3block.html" observed from the publish option of Matlab. 

We also include a real data example demonstrating the interconnection between cluster 7 and 8 of Fig. 3(f) in the paper. The code and published results are included in "laml_s3_7v8.m" and "laml_s3_7v8.html". Please use dropbox link to download the "laml_s3_step1.mat": https://www.dropbox.com/sh/iqrz1gcqp0bug5w/AAB3RCyFk28GI-KhMgVpgarba?dl=0 for the input data of "laml_s3_7v8.m".

# Implementation
Please use addpath to specify the folder of the codes. For example use: "addpath('/Users/qwu/Dropbox/Network_program-master/NICE_folder/NICE_detection')" implement the code saved under folder '/Users/qwu/Dropbox/Network_program-master/NICE_folder/NICE_detection'.
