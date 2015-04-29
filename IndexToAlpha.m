function alpha = IndexToAlpha(index,pts,alphamin,alphamax)
% Computes the regularization parameter for given index.
    
alpha = alphamin + (alphamax - alphamin)*(index-1)/(pts-1);
