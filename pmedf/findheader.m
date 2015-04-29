% function [res,pos]=findheader(header,what,kind)
% used to read header of edf-files (.edf), info-file (.info) and parameter file (.par)
% reads value of what from the string header
% what can be 'date'
% to obtain only the hour, use res(12:19)
% kind is integer, float, motor, counter or string
% origin: peter

function [res,pos]=findheader(hd,what,kind)

if nargin < 3
 	disp('3 input arguments required:')
 	help findheader
 	return
end

offset = 16;

% unused space between start of what and value of what
% contains what,'=' or ' =', and a number of blanks
% the exact offset for an edf-file is 16, for the .info file it is 24
% 16 works for both


if strcmp(what,'date')
	
	pos=findstr(hd,what);
	if ~isempty(pos)
		pos_equal=findstr(hd(pos:end),'=');
		pos_equal=pos(1)+pos_equal(1)-1;
		pos_semicolon=findstr(hd(pos_equal:end),';');
		pos_semicolon=pos_equal+pos_semicolon(1)-1;
		res=hd(pos_equal+2:pos_semicolon-2);
	else	% what is not in header
		res = [];
	end


elseif strcmp(kind,'motor')
	pos=findstr(hd,'motor_mne');
	if ~isempty(pos)
		posend=findstr(hd(pos:end),';');
		posend=pos-2+posend(1);
		pos=pos+offset+1;
		motnum=1;
		res = []; % in case motor does not exist
		while pos<posend; % looking number of motor 
			motor=sscanf(hd(pos:posend),'%s',1);
			if strcmp(motor,what) % motor was found -> taking corresponding position
				pos=findstr(hd,'motor_pos');
				pos=pos+offset+1;
				res=sscanf(hd(pos:end),'%g',motnum);
				res=res(motnum);
				break
			end
			motnum=motnum+1;
			pos=pos+1+length(motor);
		end
	else
		res = [];
	end
	
elseif strcmp(kind,'counter')
	pos=findstr(hd,'counter_mne');
	if not(isempty(pos))
		posend=findstr(hd(pos:end),';');
		posend=pos-2+posend(1);
		pos=pos+offset+1;
		motnum=1;
		res = []; % in case counter does not exist
		while pos<posend; % looking number of motor 
			motor=sscanf(hd(pos:posend),'%s',1);
			if strcmp(motor,what) % motor was found -> taking corresponding position
				pos=findstr(hd,'counter_pos');
				pos=pos+offset+1;
				res=sscanf(hd(pos:end),'%g',motnum);
				res=res(motnum);
				break
			end
			motnum=motnum+1;
			pos=pos+1+length(motor);
		end
	else
		res = [];
	end
	
else % kind is not 'motor' 'counter' 'date'
		
	pos=findstr(hd,what);
	if not(isempty(pos))
		pos = pos(1);
		pos_equal=findstr(hd(pos:end),'=');
		pos_equal=pos+pos_equal(1);

		switch kind
			case 'integer',
				res=sscanf(hd(pos_equal:end),'%i');
			case 'float',
				res=sscanf(hd(pos_equal:end),'%f');
			case 'string',
				res=sscanf(hd(pos_equal:end),'%s',1);
			case 'list',
                                pos_end = findstr(hd(pos_equal:end),'}');
				res=hd(pos_equal:pos_equal+pos_end-1);
			otherwise,
				res=[];
		end
	else
		res=[];
	end
end

  
