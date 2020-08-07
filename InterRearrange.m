function [IR,IC,C1_sort,C2_sort,C12_sort,C] = InterRearrange(C1,C2,C12,r_max,direction)
C12(abs(C12)<r_max)=0;
if direction == 'neg'
    C12_raw = C12;
    C12 = abs(C12);
end
[row,col]=size(C12);
row_orig_dens=[];
col_orig_dens=[];
for index=1:row 
    row_orig_dens(index)=sum(C12(index,:));
end

for index=1:col 
    col_orig_dens(index)=sum(C12(:,index));
end
[~,IR]=sort(row_orig_dens,'ascend');
[~,IC]=sort(col_orig_dens,'descend');
C12_sort=C12(IR,IC);
if direction == 'neg'
    C12_sort = C12_raw(IR,IC);
end
C1_sort = C1(IR,IR);
C2_sort = C2(IC,IC);
C=[C1_sort, C12_sort; C12_sort',C2_sort];
end
