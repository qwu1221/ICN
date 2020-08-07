function[dist] = KL(y1,y2,width)
x = min([min(y1) min(y2)]):width:max([max(y1) max(y2)]);
p = hist(y1,x);
q = hist(y2,x);
p = p/size(y1,2);
q = q/size(y2,2);
eps = 1e-16;
dist = sum(p.*(log2(p+eps)-log2(q+eps)));
end
