function CellTest(n,d)

fprintf('cell\n')
a=cell(1,n);
for i=1:numel(a)
    a{i}=rand(d);
end;
tic
median(cat(3,a{:}),3);
toc
clear a


fprintf('cell\n')
a=cell(1,n);
for i=1:numel(a)
    a{i}=rand(d);
end;
tic
median(cat(3,a{:}),3);
toc
clear a
