function [P,A]=KLtest(s,null,true_dist,a,width)
% s: distribution to test
% null: vector to sample null distributions
% true: distribution to compare with

% All s, null, true are row vectors
[s_kl] = KL(s,true_dist,width);
t_kl=[];
%h = waitbar(0,'Please wait...');
for i=1:1000
    %waitbar(i / 1000)
t = randsample(null,size(s,2),true);
[t_kl(i)] = KL(t,true_dist,width);
end
P = sum(t_kl>s_kl)/size(t_kl,2);
if P < a
    A=1;
else
    A=0;
end
%close(h)
end
