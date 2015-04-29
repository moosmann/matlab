if 0,
    imtool(intFT_0p01_E30_z0p01_c00000);
    imtool(intFT_0p01_E10_z0p01_c00000);
    imtool(intFT_0p01_E10_z0p02_c00000);
    imtool(intFT_0p01_E10_z0p04_c00000);
    imtool(intFT_0p01_E10_z0001_c00000);
end;
if 0,
    ishow(bro_0p01_E30_z0p01_c00000,1);
    ishow(bro_0p01_E10_z0p01_c00000,1);
    ishow(bro_0p01_E10_z0p02_c00000,1);
    ishow(bro_0p01_E10_z0p04_c00000,1);
    ishow(bro_0p01_E10_z0001_c00000,1);
end;
if 0,
    ishow(ctf_0p01_E30_z0p01_c00000,1);
    ishow(ctf_0p01_E10_z0p01_c00000,1);
    ishow(ctf_0p01_E10_z0p02_c00000,1);
    ishow(ctf_0p01_E10_z0p04_c00000,1);
    ishow(ctf_0p01_E10_z0001_c00000,1);
end;
if 1,
    ishow(abs(bro_0p01_E30_z0p01_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(bro_0p01_E10_z0p01_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(bro_0p01_E10_z0p02_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(bro_0p01_E10_z0p04_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(bro_0p01_E10_z0001_c00000(:,:,1)-pha_0p01),1);
    me(bro_0p01_E30_z0p01_c00000,pha_0p01);
    me(bro_0p01_E10_z0p01_c00000,pha_0p01);
    me(bro_0p01_E10_z0p02_c00000,pha_0p01);
    me(bro_0p01_E10_z0p04_c00000,pha_0p01);
    me(bro_0p01_E10_z0001_c00000,pha_0p01);
end;
if 1,
    ishow(abs(ctf_0p01_E30_z0p01_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(ctf_0p01_E10_z0p01_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(ctf_0p01_E10_z0p02_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(ctf_0p01_E10_z0p04_c00000(:,:,1)-pha_0p01),1);
    ishow(abs(ctf_0p01_E10_z0001_c00000(:,:,1)-pha_0p01),1);
    me(ctf_0p01_E30_z0p01_c00000,pha_0p01);
    me(ctf_0p01_E10_z0p01_c00000,pha_0p01);
    me(ctf_0p01_E10_z0p02_c00000,pha_0p01);
    me(ctf_0p01_E10_z0p04_c00000,pha_0p01);
    me(ctf_0p01_E10_z0001_c00000,pha_0p01);
end;
