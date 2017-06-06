function SET_DEFAULT( )
% Add all parameters, i.e. MATLAB variables given in the current
% workspace at the moment of calling this function to the cell struct
% 'DEFAULT'.
%
% Written by Julian Moosmann, 2017-06-02, last mod: 2017-06-03
%
% SET_DEFAULT()

evalin( 'caller', 'ADD(''default'')' )