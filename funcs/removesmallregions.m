function mask_candidates = removesmallregions(matrix_candidates,...
                                               props_candidates,...
                                               threshold,...
                                               th_cir)

%num of candidate regions
region_num = max(max(matrix_candidates));
mask_candidates = matrix_candidates;

[m,n]=size(matrix_candidates);

%total area of the image
img_area = m * n;

%if only one region found, keep it, else, pick out the regions larger than
%the threshold.

if(region_num == 1)
    return;
else
    for i = 1:region_num
        if (props_candidates(i,1).FilledArea < threshold * img_area)
           mask_candidates(mask_candidates == i) = 0;
        end
        %circularity test
        if(th_cir > 0) % 4 parameters input
             if( 4 * pi * props_candidates(i,1).FilledArea / ...
                       (props_candidates(i,1).Perimeter ^ 2) < th_cir )
                   mask_candidates(mask_candidates == i)= 0;
             end
        end
    end
end


    