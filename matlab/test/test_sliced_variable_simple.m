function test_sliced_variable_simple()

[mem_free, mem_avail_cpu, mem_total_cpu] = free_memory;

fprintf( '\n Create array' )
a = ones( [2001,2002, 2000]);

fprintf( '\n Call myfunc' )
a = myfunc(a);

fprintf( '\n End' )
end

function a = myfunc(a)
   fprintf( '\n Start loop' )
   s1 = size(a,1);
   s2 = size(a,2);
    for n = 1:size(a,3)
        % Do something
        im = a(:,:,n);
        im = im + 1;                
        a(:,:,n) = im(1:s1,1:s2);
        if n ==1
            fprintf( '\n First iteration finished' )
        end 
    end
end