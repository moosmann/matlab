function index = AlphaToIndex(alpha,alphamin,alphamax,pts)
    
index = (alpha - alphamin)/(alphamax - alphamin)*(pts-1) + 1;
