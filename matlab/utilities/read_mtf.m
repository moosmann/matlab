function read_mtf(path)
% Display mtf image

if nargin < 1
    path = pwd;
end

ca

%% Read images
fn = sprintf('%s/*.da0',path);
d = dir(fn);
fn = sprintf('%s/%s',d.folder,d.name);
dar = read_dat_jm(fn);
dar = FilterPixel(dar);
dar = Binning(dar,2);

fn = sprintf('%s/*.re0',path);
d = dir(fn);
fn = sprintf('%s/%s',d.folder,d.name);
ref = read_dat_jm(fn);
ref = FilterPixel(ref);
ref = Binning(ref,2);

fn = sprintf('%s/*.ob0',path);
d = dir(fn);
fn = sprintf('%s/%s',d.folder,d.name);
obj = read_dat_jm(fn);
obj = FilterPixel(obj);
obj = Binning(obj,2);

im = (obj - dar) ./ (ref -dar);
ym = mean(im,2);
s = floor(size(im,2)/2);
y1 = im(:,s);

ishow(dar)
%xticks('auto'),yticks('auto')
ishow(ref)
ishow(obj)
ishow(im)

%dmax = max(dar(:));
figure('Name','MTF')
x = 1:numel(ym);

plot(x,ym,x,y1)
legend({'mean',sprintf('y = %u',s)})
axis tight

fprintf('\n')