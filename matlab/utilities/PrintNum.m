function PrintNum(nn,newLine,printFormat)
% Print input number 'nn' with Format 'printFormat' to Command Window
% output and print a line when mod(nn,newLine).
%
%Written by Julian Moosmann,2013-10-09

%% Default parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
    newLine = 30;
end
if nargin < 3
    printFormat = '%4u ';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(printFormat,nn); 
if mod(nn,newLine)==0
    fprintf('\n   ');
end