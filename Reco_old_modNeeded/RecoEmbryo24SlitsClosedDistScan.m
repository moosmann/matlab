function RecoEmbryo24SlitsClosedDistScan(intensity,dist,ComputeCorrection)

if nargin<3, ComputeCorrection=0;end

Edp = [25 1 1e-6];

int=intensity(:,:,dist);
phistack=int;
phicorstack=int;
alpha = cell(11,1);
%alpha{1}=[0.5 1 1.5 2 2.5 3.0 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0 4.5 5.5 6 6.5];
if isempty(alpha{dist}),
    if ComputeCorrection==0
        alpha=[0.5:0.25:5.5];
    else
        alpha=[0.5:0.5:6];
    end
else
    alpha=alpha{dist};
end;
N=length(alpha);
for i=[1:N],
    phi = Reco(int,alpha(i),Edp,ComputeCorrection);
    if ComputeCorrection==0,
        phistack=cat(3,phistack,phi(:,:,1));
    else
        phistack=cat(3,phistack,phi(:,:,1));
        phicorstack=cat(3,phicorstack,phi(:,:,2));
    end
end
if ComputeCorrection==0
    assignin('base',['phi' num2str(dist,'%02u')],phistack);
else
    assignin('base',['phi' num2str(dist,'%02u')],phistack);
    assignin('base',['phicor' num2str(dist,'%02u')],phicorstack);
end
