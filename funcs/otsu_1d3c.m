%segment img_input into 3 classes based on the 1D histogram
%a simple extension of ordinary 2 class 1D Otsu method
%
%kuanli@cityu.edu.hk, GPL V2, 2011.6
%@input: img_input (uint8 image)
%@output:
%(1)img_res
%(2)th1, threshold 1,the smaller one
%(3)th2, threshold 2,the bigger one


function [img_res,th1,th2] = otsu_1d3d(img_input)

format long

[m,n] = size(img_input);

hist = zeros(1,256);%histogram

%notice: [0,255] -> [1,256]
for i = 1:m
    for j = 1:n
        hist(img_input(i,j)+1) = hist(img_input(i,j)+1) + 1;
    end
end

p = zeros(1,256);
total_avg_intensity = 0;

for i = 1:256
    p(i) = hist(i)/(m*n);
    total_avg_intensity = total_avg_intensity + i * p(i);
end

%-----------------------------------------
err = 0;
errmax = 0;

th1 = 0;
th2 = 0;

w1=0;
w2=0;
w3=0;

tmp1=0;
tmp2=0;
tmp3=0;

u1=0;
u2=0;
u3=0;

for t1 = 1:256
    if(hist(t1)==0)
        continue;
    end
    
    w1 = w1+p(t1);
    tmp1 = tmp1 + t1*p(t1);
    u1 = tmp1/w1;
    
    w2=0;
    tmp2=0;
    
    w3=0;
    tmp3=0;
    
    for t2 = (t1+1):256
        tmp3 = tmp3 + t2*p(t2);
    end
    
    for t2 = (t1+1):256
        if(hist(t2)==0)
            continue;
        end
        
        w2=w2+p(t2);
        tmp2 = tmp2+t2*p(t2);
        u2 = tmp2/w2;
        
        w3 = 1 - w1 - w2;
        tmp3 = tmp3 - t2*p(t2);
        
        if(w3 == 0)
            continue;
        end
        
        u3 = tmp3/w3;
        
        err = w1*(u1-total_avg_intensity)^2+...
              w2*(u2-total_avg_intensity)^2+...
              w3*(u3-total_avg_intensity)^2;
        
        if(err > errmax)
            errmax = err;
            th1 = t1;
            th2 = t2;
        end
    end
end

th1 = th1 -1;
th2 = th2 -1;


img_res = zeros(m,n);

for i = 1:m
    for j = 1:n
        
        if(img_input(i,j) >= th2)
                %img_res(i,j) = th2;
                img_res(i,j) = 200;
        else if(img_input(i,j) > th1 )
                %img_res(i,j) = th1;
                img_res(i,j) = 100;
            end
        end
    end
end

img_res = uint8(img_res);
