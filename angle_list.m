
N_outer = 701;
N_fovs = 15;

fov_angle_lists = (0:(N_outer-1)).*ones(N_fovs,1);

for k = 0:N_outer-1
    keepfactor = 8-abs(k-7);
    fov_angle_lists(k,1:keepfactor:end)=0;
end


counter = 1;
for ii = 1:N_outer
    if fovfov_angle_lists(k,8) == 1
        for k = 1:N_fovs
            final_list(counter) = fov_angle_list(k,ii);
            fov_angle_list(k,ii) = 0;
            counter = counter+1;
        end  
    end
end

final_list = [final_list,fov_angle_list(ne(fov_angle_list,0))];


N_inner = 100;

inner_list = (360/N_inner)*(0:N_inner-1);
counter = 1;
counter2 = 1;
for k = 0:N_fovs-1
    keepfactor = 8-abs(k-7);
    this_list = (360/(N_inner*keepfactor))*(0:(N_inner*keepfactor)-1);
    for ii = 1:length(this_list)
        if any(this_list(ii) == inner_list)
            final_list(counter) = this_list(ii);
            counter = counter + 1;
        else
            save_for_later(counter2) = this_list(ii);
            counter2 = counter2 +1;
        end
    end
end



%%

N_fov = 15;
N_middle = 1000;

angles_middle = (360/N_middle)*(0:N_middle-1);

for k = 0:N_fov-1
    multiplier = 8-abs(k-7);
    this_list = (360/(N_inner*multiplier))*(0:(N_inner*multiplier)-1);

end
%% angles not in the first scan

N_fov = 15;
N_middle = 1000;

angles_middle = (360/N_middle)*(0:N_middle-1);

counter = 1;
counter2 = 1;
for k = 0:N_fovs-1
    multiplier = abs(k-7);
    this_list = (360/(N_inner*multiplier))*(0:(N_inner*multiplier)-1);
    for ii = 1:length(this_list)
        if any(this_list(ii) == angles_middle)
            first_list(counter) = this_list(ii);
            counter = counter + 1;
        else
            second_list(counter2) = this_list(ii);
            counter2 = counter2 +1;
        end
    end
end




















