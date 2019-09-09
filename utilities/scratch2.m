clear all

assign_default( 'a.b', 1 )
disp( a.b )
assign_default( 'a.b', 2 )
disp( a.b )

assign_default( 'b', 3 )
disp( b )
assign_default( 'b', 4 )
disp( b )

assign_default( 'c', [] )
disp( c )
assign_default( 'c',5 )
disp( c )

assign_default( 'c',5, 1 )
disp( c )

