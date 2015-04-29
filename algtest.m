readdat;

alpha = 12;
phn = Reco(dat,alpha,1);
pho = rec(dat,1,alpha);

phn1 = phn(:,:,1);
phn2 = phn(:,:,2);
pho1 = pho(:,:,1);
pho2 = pho(:,:,2)+pho(:,:,3)+pho(:,:,4);

domain(phn1);domain(pho1);
domain(phn2);domain(pho2);

