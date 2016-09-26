for rotax=[RotAxisPos 985.4 985.6 986.1 987 1007 1007.1 1007.6 1008 1008.1 1008.9 1023 1023.1 1023.9 1024 1024.1 1024.9 1025]
    imPad=RotAxisSymmetricPadding(im,rotax);
    im2Pad=RotAxisSymmetricPadding(im2,rotax);
    out=ImageCorrelation(imPad,fliplr(im2Pad),0,0);
    fprintf('Rotation axis: %6.1f, Output images: dimHor: %g %g, RotAxis: %f\n',...
        rotax,size(imPad,2),size(im2Pad,2),out.VerticalRotationAxisPosition);
end