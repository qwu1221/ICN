function [RR Clu_size fp fn fp_intra fn_intra fp_inter fn_inter] = sim_edge2(Cor_est,Cor_true,c,true_cut,true_clu)
%% This function is to calculate edge FPR and FNR for intra- and inter-community separately
%  Cor_est:        a sample correlation matrix
%  Cor_true:   a true covariance matrix/correlation matrix
%  c:              number of communities for true matrix (only take care the first c detected cluster for Cor_est)
%  true_cut:       cut-off for true significant edges from true correlation
%  matrix


%% Get true significant edges from true correlation matrix by cutting off
n=size(Cor_true,1);
if sum(sum(diag(Cor_true)))>0
Cor_true = Cor_true-eye(n);
end
true_edges = (abs(Cor_true)>true_cut);
%figure;imagesc(true_edges);colormap jet;colorbar;snapnow
true_intra = zeros(n);
for i=1:c 
true_intra(true_clu{i},true_clu{i}) = ones(size(true_clu{i},2))-eye(size(true_clu{i},2));
end
%figure;imagesc(true_intra);colormap jet;colorbar;snapnow
true_inter = (true_edges-true_intra)>0;
%figure;imagesc(true_inter);colormap jet;colorbar;snapnow


%% Perform permutation on sample correlation matrix
addpath('/Users/qwu/Downloads/Don/Interconnected')
addpath('/Users/qwu/Downloads/Network_program-master/NICE_folder/NICE_detection')
for i=1:n
    Cor_est(i,i)=0;
end
perm_matrix = squareform(1:(n*(n-1)/2));
node_perm_idx = randperm(n);
perm_matrix = perm_matrix(node_perm_idx,node_perm_idx);
perm_vec = squareform(perm_matrix);

Cor_vec = squareform(Cor_est);
Cor_vec = Cor_vec(perm_vec)';
Cor_perm = squareform(Cor_vec);

%% NICE detection
[CindxVICC,CIDVICC,ClistVICC]=NICE(Cor_vec, 0.7, 0, 10);
Cor_sort = Cor_perm(ClistVICC,ClistVICC);

%% Find non-singleton clusters (actually size >5)
Clu_size = [];
for i=1:size(CIDVICC,2)
    Clu_size(i) = size(find(CindxVICC==CIDVICC(i)),2);
end

[large_clu,large_clu_idx]=sort(Clu_size,'descend');
Clu_size = Clu_size(Clu_size>2); k = size(Clu_size,2);
% If the number of non-singleton clusters are smaller than c, we set this
% detection as 'error'
clear RR
if k==1
    RR = 'error';
    fp_intra = 0;
    fn_intra = 1;
    fp_inter = 0;
    fn_inter = 1;
    fp = 0;
    fn = 1;
    return
end
if k>1 & k<c
     c=min(k,c);
end
addpath('/Users/qwu/Downloads/Don/Interconnected')
for i=1:k
A{i} = find(CindxVICC==CIDVICC(large_clu_idx(i)));
end
idx = 1:size(Cor_perm);
idx_select = [A{1:k}];
idx_left = setdiff(idx,idx_select);

clu00 = Cor_perm(idx_left,idx_left); %% random-graph as singleton part

%% Diagnal blocks (largest c)
for i=1:c
Diag{i}=Cor_perm(A{i},A{i});
end

%% Off-diagnals blocks
clear Off Off_vec
Off_1vec = [];
for i=1:(c-1)
    for j=(i+1):c
        Off{i,j} = Cor_perm(A{i},A{j});
        CC = Cor_perm(A{i},A{j});
        Off_vec{i,j} = CC(:);
        VV = CC(:);
        Off_1vec = [Off_1vec VV'];
    end
end


%% Test interconnectivity by KL
Off_2 = Cor_perm([A{1:k}],idx_left);
Off_2vec = Off_2(:);
true_dist = squareform(clu00);
null = [Off_1vec Off_2vec' true_dist];
width= 0.001;

for i=1:(c-1)
    for j=(i+1):c
        s = Off_vec{i,j}';
        %histogram(s,'Normalization','probability'); hold on;histogram(null,'Normalization','probability'); hold on; histogram(true_dist,'Normalization','probability')
        [P,R]=KLtest(s,null,true_dist,0.01,width);
        %[i j P R]
        RR(i,j) = R;
    end
end
%RR

%% Get significant interconnected edges by cutting-off from 'InterCut'
clear raw
B= zeros(n);
for i=1:c 
    raw{i} = node_perm_idx(A{i});
B(raw{i},raw{i})=ones(size(A{i},2))-eye(size(A{i},2));
end
%figure;imagesc(true_intra);colormap jet;colorbar;snapnow
est_intra_vec = squareform(B);
true_intra_vec = squareform(true_intra);

fp_intra = sum(sum(est_intra_vec>true_intra_vec))/(size(true_intra_vec,2)-sum(sum(true_intra_vec)));
fn_intra = sum(sum(est_intra_vec<true_intra_vec))/sum(sum(true_intra_vec));
%fp_intra = sum(sum(est_intra_vec>true_intra_vec));
%fn_intra = sum(sum(est_intra_vec<true_intra_vec));


C = zeros(n);

if sum(sum(RR))>0
[s,t,u]=find(RR);

clear C_cut
r_cut = [];
lambda0=0.6;
r=0.1:0.005:0.8;
for i=1:size(s,1)
       C1 = Diag{s(i)};
       C2 = Diag{t(i)};
       C12 = Off{s(i),t(i)};
       r_max=InterCut(C1,C2,C12,r,lambda0);
       r_cut(i) = r_max(1);
       C_cut{i} = abs(C12)>r_cut(i);
       C(raw{s(i)},raw{t(i)})=abs(C12)>r_cut(i);
       C(raw{t(i)},raw{s(i)})=abs(C12')>r_cut(i);
end 
end

%figure;imagesc(C);colormap jet;colorbar;snapnow
%figure;imagesc(true_inter);colormap jet;colorbar;snapnow

est_inter_vec = abs(squareform(C));
true_inter_vec = squareform(true_inter);

%% Compare true and sample by significant edges
fp_inter = sum(sum(est_inter_vec>true_inter_vec))/(size(true_inter_vec,2)-sum(sum(true_inter_vec)));
fn_inter = sum(sum(est_inter_vec<true_inter_vec))/sum(sum(true_inter_vec));
%fp_inter = sum(sum(est_inter_vec>true_inter_vec));
%fn_inter = sum(sum(est_inter_vec<true_inter_vec));
est_vec = squareform(B+C);
true_vec = squareform(true_edges);

%% Compare true and sample by significant edges
fp = sum(sum(est_vec>true_vec))/(size(true_vec,2)-sum(sum(true_vec)))
fn = sum(sum(est_vec<true_vec))/sum(sum(true_vec))

end
