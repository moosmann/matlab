function VarOut = VarDefault2(VarName,VarDefaultValue)
% Check if variable 'VarName' exists. If not assign default value
% 'VarDefaulValue' to it'
%
% VarName: String
% VarDefaultValue: any


if exist(VarName,'var')
   VarOut = 0;
else
    fprintf('else')
    VarOut = VarDefaultValue;
end

