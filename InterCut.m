function [r_max]=InterCut(C1,C2,C12,r,lambda0)
obj=[];
for i=1:size(r,2)
W=C12(abs(C12)>r(i));
num = sum(abs(W))+sum(abs(squareform(C1)))+sum(abs(squareform(C2)));
denom = size(W,1)+size(squareform(C1),2)+size(squareform(C2),2);
obj(i) = num/(denom)^lambda0;
end
r_max = r(obj==max(obj));
end