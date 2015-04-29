#------------------
# creates the name of a single reference file
# no device server: ref10000.edf ...
# with device server: ref0000_0000.edf ...

function y=refsname(number,pos,devserver)
    if devserver
        y = sprintf('ref%4.4i%s%4.4i.edf',number-1,'_',pos);
    else
        y = sprintf('ref%1i%4.4i.edf',number,pos);
    end
endfunction
