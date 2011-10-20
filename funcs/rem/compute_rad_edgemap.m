function rem = compute_rad_edgemap(ori,centroid, Options)

assert(isa(ori,'double'),...
 'image input when computing the rad edge map is supposed to be double type');

[m,n]=size(ori);

interpretY = [];
interpretX = [];
interpretValue = [];

rem = zeros(m,n);

cen(1) = round(centroid(2));
cen(2) = round(centroid(1));

%boundary pixels
tmp = ones(m,n);
boundary = bwboundaries(tmp);
boundary_pixels = boundary{1,1};

wait_bar = waitbar(0,'computing radiating edge map...');

pixel_num = 2 *(m+n) - 4 ;

for i = 1 : pixel_num 
    
    if ~isempty(wait_bar)
        waitbar(i/pixel_num,wait_bar);
    end

    
    [value,pos]=get_single_line(ori,cen,[boundary_pixels(i,1),...
                                         boundary_pixels(i,2)]);
    forwardDiff = compute_forward_difference(value);
    
    if Options.isstack_refine
        segInfo = compute_segment(forwardDiff, Options.stackrefine_theta);
        forwardDiff = stack_refine(forwardDiff, segInfo);
    end
    
    if Options.ispositive_suppress
        forwardDiff = positive_suppress(forwardDiff, ...
                                       Options.positivesuppress_ratio);
    end
    
    radiating_gradient = compute_gradient(forwardDiff);
   
    %Normalize
    minRG = min(radiating_gradient);
    maxRG = max(radiating_gradient);
    
    normRG = 255*(radiating_gradient - minRG)/(maxRG - minRG);

    interpretY = [interpretY; pos(:,1)];
    interpretX = [interpretX; pos(:,2)];
    interpretValue = [interpretValue; normRG];
    
end


F = TriScatteredInterp(interpretX,interpretY,interpretValue);

for i = 1 : m
    for j = 1 : n
        rem(i,j) = F(j,i);
        if(isnan(rem(i,j)))
            rem(i,j) = 0;
        end
    end
end

close(wait_bar);






