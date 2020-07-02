% function str=headerstring(what,value,kind)
% used by writeheader
% returns a string for incorporation in an edf-header
% what is the information item (e.g. 'Dim_1')
% value of the item
% kind is string, integer or float
% origin: peter

function str=headerstring(what,value,kind)

offset=16;
if length(what)<offset
	str=[what blanks(offset-1-length(what)) '= '];
	switch kind
		case 'string',
		   	str=[str sprintf('%s ;\n',value')];
		case 'integer',
         		str=[str sprintf('%i ;\n',value')];
      		case 'float'
         		str=[str sprintf('%g ;\n',value')];
	end
end
