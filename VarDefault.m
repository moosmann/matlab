function VarDefault(VarName,VarDefaultValue)
% Check if variable 'VarName' exists. If not assign default value
% 'VarDefaulValue' to it'
%
% VarName: String
% VarDefaultValue: any

if ~exist(VarName,'var')
    assignin('caller',VarName,VarDefaultValue)
end



