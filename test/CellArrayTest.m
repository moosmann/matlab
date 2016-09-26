function CellArrayTest()

clear all

n=10;
d=500;

fprintf(1,'\nstart comparison\n')

fprintf(1,'cell\n')
a=cell(1,n);
for i=1:numel(a)
    a{i}=rand(d);
end;
tic
median(cat(3,a{:}),3);
toc
clear a


fprintf(1,'array\n')
b=rand(d,d,n);
tic
median(b,3);
toc
clear b


fprintf(1,'cell\n')
a=cell(1,n);
for i=1:numel(a)
    a{i}=rand(d);
end;
tic
median(cat(3,a{:}),3);
toc
clear a


fprintf(1,'array\n')
b=rand(d,d,n);
tic
median(b,3);
toc
clear b


fprintf(1,'cell\n')
a=cell(1,n);
for i=1:numel(a)
    a{i}=rand(d);
end;
tic
median(cat(3,a{:}),3);
toc
clear a



clear all