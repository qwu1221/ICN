# ICN
Extracting Interconnected Communities in Gene co-expression Networks

Step 1: detect communities described in step 1 of manuscript using NICE.m;

Step 2: test the interconnected community pairs using algoirthm 1 of step 2 by KLtest.m;

Step 3: identify connecting edges as algorithm 2 of step 3 by InterCut.m;

sim_edge2.m calculates the edge-wise FPR and FNR with a given true correlation matrix and observed correlation matrix;

demo_3block.m is a small simulation example with ICN detection;

laml_s3_7v8.m is a the real data example demonstrating the interconnection between cluster 7 and 8. Please use dropbox link to download the 'laml_s3_step1.mat': https://www.dropbox.com/sh/iqrz1gcqp0bug5w/AAB3RCyFk28GI-KhMgVpgarba?dl=0
