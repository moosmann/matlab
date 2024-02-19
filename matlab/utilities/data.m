p = '/asap3/petra3/gpfs/p05/2021/data';
p = '/asap3/petra3/gpfs/p07/2023/data';
p = '/asap3/petra3/gpfs/p05/2023/data';

pp = [p '/*/processed/*/reco*/*/float_rawBin*'];
pa = [p '/*/processed/*/reco/float_rawBin*'];
pd = [p '/*/processed/*/Reco/reco_phase/Int*'];

strp = sprintf('ls -d1 %s',pp);
stra = sprintf('ls -d1 %s',pa);
strd = sprintf('ls -d1 %s',pd);

fprintf('\n%s\n',strp)
fprintf('\n%s\n',stra)
fprintf('\n%s\n',strd)

[~,rd] = unix(strd);
cd = textscan(rd,'%s');
cd = cd{1};
for n = numel(cd):-1:1
    c = cd{n};
    [fp,name,ext] = fileparts(c);
    %fprintf('\n %s %s %s',fp,name,ext)
    if ~isempty(ext)
        cd(n) = [];
 %       fprintf('\n %s %s %s',fp,name,ext)
    end
    if contains(name,'recoBin')
        cd(n) = [];
%        fprintf('\n %s %s %s',fp,name,ext)
    end
end
c = cd;
mtotal = 0;
for n = 1:numel(c)
    if n == 2701
        continue
    end
    s = sprintf('du -hs %s',c{n});
    [~,r] = unix(s);
    %fprintf(r)
    t = textscan(r,'%fG %s');
    m = t{1};
    mtotal = mtotal + m;
    fprintf('%6u %10.0f %s',n,mtotal,r)
end








[~,rp] = unix(strp);
[~,ra] = unix(stra);

cp = textscan(rp,'%s');
cp = cp{1};

ca = textscan(ra,'%s');
ca = ca{1};



whos cp ca cd

for n = numel(cp):-1:1
    c = cp{n};
    [fp,name,ext] = fileparts(c);
    %fprintf('\n %s %s %s',fp,name,ext)
    if ~isempty(ext)
        cp(n) = [];
%        fprintf('\n %s %s %s',fp,name,ext)
    elseif contains(name,'recoBin')
        cp(n) = [];
 %       fprintf('\n %s %s %s',fp,name,ext)
    end
end

for n = numel(ca):-1:1
    c = ca{n};
    [fp,name,ext] = fileparts(c);
    %fprintf('\n %s %s %s',fp,name,ext)
    if ~isempty(ext)
        ca(n) = [];
 %       fprintf('\n %s %s %s',fp,name,ext)
    end
    if contains(name,'recoBin')
        ca(n) = [];
%        fprintf('\n %s %s %s',fp,name,ext)
    end
end

c = cat(1,ca,cp,cd);

whos cp ca c

mtotal = 0;
for n = 1:numel(c)
    if n == 2701
        continue
    end
    s = sprintf('du -hs %s',c{n});
    [~,r] = unix(s);
    %fprintf(r)
    t = textscan(r,'%fG %s');
    m = t{1};
    mtotal = mtotal + m;
    fprintf('%10.0f %s',mtotal,r)
end

