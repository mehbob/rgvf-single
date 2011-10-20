function [res, props] = getmask(type,ori,ths,img,ratio, th_cir)

    if nargin == 5
        th_cir = 0;
    end
    
    assert(isa(img,'uint8'),...
            'img is supposed to be uint8 type');

    %close operation when type = 'cyto'
    if strcmp(type,'cyto')
        imgBin = im2bw(ori,ths(2)/255.0);
        imgBin = bwfill(imgBin, 'holes', 8);
        %closeStructure = strel('disk',1);
        %imgBin = imclose(imgBin,closeStructure);
    else if strcmp(type,'nu')
            imgBin = im2bw(ori,ths(3)/255.0);
            imgBin = bwfill(imgBin, 'holes', 8);
         else
             error('type should be either "cyto" or "nu"');
         end
    end
    
    %label the binary image
    mat =bwlabel(imgBin,8);

    %get the properties of the regions
    props  = regionprops(mat, img, 'all');
    
    res = removesmallregions(mat, props, ratio, th_cir);
    