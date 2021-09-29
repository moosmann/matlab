function test_sliced_variable_simple( vol_mem )

if nargin < 1
    vol_mem = 10; % GB;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ca
[mem_free, mem_avail_cpu, mem_total_cpu] = free_memory;
mf0 = mem_free;
ma0 = mem_avail_cpu;
[mm,p] = getmem;
ms = getrss(p);
fprintf( '\n RAM:' )
fprintf( '\n free      : %.0f GiB (%.1f%%)', mem_free/1024^3, 100 * mem_free/mem_total_cpu )
fprintf( '\n available : %.0f GiB (%.1f%%)', mem_avail_cpu/1024^3, 100*mem_avail_cpu/mem_total_cpu )
fprintf( '\n total     : %.0f GiB', mem_total_cpu/1024^3)
fprintf( '\n MATLAB : %f (%f) GiB', mm / 100 * mem_total_cpu / 1024^3, ms / 1024^3)

s = round( (vol_mem * 1024^3 / 4 )^(1/3) );
a_size = s + [0 ,1, 2];
a_mem = prod( a_size ) * 4;
fprintf( ' \n\n volume size : %u, %u, %u', a_size )
fprintf( ' \n volume memory : %.0f B = %.2f GiB ', a_mem, a_mem / 1024^3)

fprintf( '\n Create array:' )
tic;
a = ones(  a_size, 'single');
t = toc;
fprintf( ' %f s = %f min', t, t/60 )

a_mem_is = Bytes( a);
fprintf( ' \n volume memory : %.0f B = %.2f GiB ', a_mem_is, a_mem_is / 1024^3)

[mem_free, mem_avail_cpu, mem_total_cpu] = free_memory;
fprintf( '\n\n RAM:' )
fprintf( '\n free      : %.0f GiB (%.1f%%)', mem_free/1024^3, 100 * mem_free/mem_total_cpu )
fprintf( '\n available : %.0f GiB (%.1f%%)', mem_avail_cpu/1024^3, 100*mem_avail_cpu/mem_total_cpu )
fprintf( '\n MATLAB : %f (%f) GiB', getmem(p) / 100 * mem_total_cpu / 1024^3, getrss(p) / 1024^3)

fprintf( '\n\n Call myfunc' )
t = toc;
[a, mm, ms, mf, ma] = myfunc(a, p);
fprintf( '\n myfunc finished: %.0f s = %.1f min',  toc -t, (toc - t) / 60)

figure('Name', 'Memory consumption')
p0 = mm';
p1 = ms';
p2 = flipud(( ma0 - sort(ma))');
p3 = flipud(( mf0 - sort(mf))');
plot( [p0, p1, p2, p3]/1024^3)
axis tight
xlabel( 'iteration #')
ylabel( 'memory / GiB')
legend( {'ps %mem', '%ps rss',  'free', 'avail'} )

clear all
[mem_free, mem_avail_cpu, mem_total_cpu] = free_memory;
fprintf( '\n\n RAM:' )
fprintf( '\n free      : %.0f GiB (%.1f%%)', mem_free/1024^3, 100 * mem_free/mem_total_cpu )
fprintf( '\n available : %.0f GiB (%.1f%%)', mem_avail_cpu/1024^3, 100*mem_avail_cpu/mem_total_cpu )

fprintf( '\n')
end

function [a, mm, ms, mf, ma] = myfunc(a, p)
   fprintf( '\n Start loop' )
   s1 = size(a,1);
   s2 = size(a,2);
   np = round(size(a,3) / 4 );
   mf = zeros( [1,size(a,3)] );
   ma = mf;
   mm = mf;
   ms = mf;
    for n = 1:size(a,3)
        % Do something
        im = a(:,:,n);
        im = im + 1;                
        a(:,:,n) = im;
        [mem_free, mem_avail_cpu, mem_total_cpu] = free_memory;            
        mf(n) = mem_free;
        ma(n) = mem_avail_cpu;
        mm(n) = getmem(p) / 100 * mem_total_cpu;
        ms(n) = getrss(p);
        %if mod(n,np) == 1         
         %   fprintf( '\n Iteration %4u: %.0f GiB, %.0f GiB, %0.f GiB', n, mem_free/1014^3, mem_avail_cpu/1024^3,  mm(n)/1024^3)
        %end 
    end
end


