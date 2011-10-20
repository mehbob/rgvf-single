function img_after_kmeans = kmeans_spatial(img_ori)
%divide the original image(@para: img_ori) into three classes
%using kmeans with spatial information

%copyright @ GPL v2, kuanli

[m,n]=size(img_ori);

img_ori_double = double(img_ori);

%para, neighborhood size (2*h+1)*(2*h+1)
h=2;

img_mean =zeros(m,n);
img_median = zeros(m,n);


%cal the mean and median value of the neighborhood
for i=1:m
    for j=1:n
        %tmp array, help to find the median value
        array_tmp = zeros(1,((2*h+1)^2));
        index_tmp = 1;
        for k=-h:h
            for w=-h:h;
                p=i+k;
                q=j+w;
                if (p<=0)|( p>m)
                    p=i;
                end
                if (q<=0)|(q>n)
                    q=j;
                end
                img_mean(i,j)=img_mean(i,j)+img_ori_double(p,q)/((2*h+1)^2);
                array_tmp(index_tmp) = img_ori_double(p,q);
                index_tmp = index_tmp + 1;
            end
        end
    img_median(i,j) = median(array_tmp);
    end
end


array_img_ori = reshape(img_ori_double,m*n,1);
img_spatial_features = [array_img_ori,...
                        reshape(img_median, m*n, 1),...
                        reshape(img_mean, m*n,1)];


%find the intial three centroids
array_unique = unique(array_img_ori);

centroid_min = min(array_unique);
centroid_max = max(array_unique);

elements_left_index = array_unique > centroid_min &...
                     array_unique < centroid_max;
                 
elements_left = array_unique(elements_left_index);
centroid_median = median(elements_left);

%call matlab func kmeans
[array_id,centroids] = kmeans(img_spatial_features,3,...
                        'Start',[centroid_min,centroid_min,centroid_min;...
                                 centroid_median,centroid_median,centroid_median;...
                                 centroid_max,centroid_max,centroid_max]);
matrix_id = reshape(array_id,m,n);

matrix_after_kmeans = zeros(m,n);

[~, id] = sort(centroids(:,1));
for i = 1:m
    for j = 1:n
        if(matrix_id(i,j)==id(1))
            %matrix_after_kmeans(i,j)=centroids(id(1),1);
            matrix_after_kmeans(i,j)=0;
        else if(matrix_id(i,j) == id(2))
                %matrix_after_kmeans(i,j) = centroids(id(2),1);
                matrix_after_kmeans(i,j) = 100;
            else 
                %matrix_after_kmeans(i,j) = centroids(id(3),1);
                matrix_after_kmeans(i,j) = 200;
            end
        end
    end
end

img_after_kmeans = uint8(matrix_after_kmeans);


