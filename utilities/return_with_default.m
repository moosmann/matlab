function variable =  return_with_default( variable, default )

if isempty( variable )
    variable = default;
end