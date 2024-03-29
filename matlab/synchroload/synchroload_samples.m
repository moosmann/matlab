function out = synchroload_samples( )
%% Ti
% CPD
ti_c = {
'1L'
'1R'
'4R'
'7L'
'22L'
'34R'
'39R'
'42L'
'9L'
'10L'
'11L'
'12L'
'14L'
'15L'
'16L'
'19L'
'25L'
'26L'
'27L'
'28L'
'29L'
'30L'
'31R'
'38R'
'41L'
'43L'
'45L'
'46L'
'48R'
'51R'
'52R'
};
ti_c = sort( ti_c );
% LOAD
ti_l = {
'13L'
'17L'
'20L'
'21L'
'23L'
'24L'
'32L'
'35R'
'37R'
'3R'
'40R'
'47R'
'49R'
'54L'
'5R'
'61R'
'6R'
'8R'
};
ti_l = sort( ti_l );
out.Ti = cat( 1,  ti_c, ti_l );
out.cpd.Ti = ti_c;
out.load.Ti = ti_l;
%% PEEK
% CPD
peek_c = {
'39L'
'21R'
'52L'
'38L'
'10R'
'28R'
'48L'
'3L'
'4L'
'8L'
'9R'
'12R'
'13R'
'15R'
'16R'
'20R'
'25R'
'26R'
'27R'
'29R'
'30R'
'35L'
'37L'
'41R'
'43R'
'45R'
'46R'
'61L'
};
peek_c = sort( peek_c );
%LOAD
peek_l = {
'11R'
'17R'
'18R'
'22R'
'23R'
'2L'
'31L'
'32R'
'36L'
'40L'
'42R'
'44R'
'47L'
'49L'
'50L'
'54R'
'5L'
'7R'
};
peek_l = sort( peek_l );
out.cpd.PEEK = peek_c;
out.load.PEEK = peek_l;
out.PEEK = cat( 1, peek_c, peek_l );
%% Mg5Gd
mg5gd_c = {
'89L'
'101L'
'78L'
'66L'
'75L'
'63L'
'69L'
'80L'
'100L'
'59R'
'62L'
'64L'
'70L'
'72L'
'85R'
'87R'
'88R'
'90R'
'91R'
'99L'
'97L'
'74L'
'76L'
'77L'
'82L'
'93R'
'94R'
'95R'
'96R'
};
mg5gd_c = sort( mg5gd_c );
% LOAD
mg5gd_l = {
'102L'
'103L'
'57R'
'60R'
'67R'
'73L'
'83L'
'84R'
'86R'
'92R'
};
mg5gd_l = sort( mg5gd_l );
out.cpd.Mg5Gd = mg5gd_c;
out.load.Mg5Gd = mg5gd_l;
out.Mg5Gd = cat( 1, mg5gd_c, mg5gd_l );

%% Mg10Gd
% CPD
mg10gd_c = {
'80R'
'74R'
'69R'
'62R'
'88L'
'87L'
'53L'
'66R'
'82R'
'99R'
'64R'
'70R'
'56L'
'60L'
'68R'
'90L'
'100R'
'97R'
'91L'
'101R'
'73R'
'77R'
'78R'
'79R'
'93L'
'95L'
'96L'
'100AL'
};
mg10gd_c = sort( mg10gd_c );
% LOAD
mg10gd_l = {
'102R'
'103R'
'57L'
'63R'
'67L'
'71R'
'72R'
'81R'
'83R'
'85L'
'89R'
};
mg10gd_l = sort( mg10gd_l );
out.cpd.Mg10Gd = mg10gd_c;
out.load.Mg10Gd = mg10gd_l;
out.Mg10Gd = cat( 1, mg10gd_c, mg10gd_l );
