ca
bottom = min([intEmean(:); intMean(:)]);
top = max([intEmean(:); intMean(:)]);

h = figure;
ax(1) = subplot(1,2,1);
h1 = imagesc(intEmean);
colormap(h,gray);
caxis([bottom top])
ax(2) = subplot(1,2,2); 
h2 = imagesc(intMean);
colormap(h,gray);
caxis([bottom top])
axis(ax,'off');
axis(ax,'image');
colorbar('peer',ax(2),[0.965064102564103 0.2 0.032051282051282 0.632]);
