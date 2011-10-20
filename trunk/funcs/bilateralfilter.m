function FA = bilateralfilter(A, size, std_c, std_s)
%BILATERALFILTER Filters a gray level image with a bilateral filter.
%   FA = BILATERALFILTER(A, size, std_c, std_s) filters the
%   gray level image A using a window of [size(1) size(2)] with a
%   standard bilateral filter. 
%   Bilateral filtering smooths images while preserving edges in contrast.
%   Therefor an adaptive filtering kernel for each image element in
%   the convolution is created. The kernel is the product of a gaussian 
%   kernel (closeness function) and a gaussian weighted similarity function
%   for pixel intensities. std_c and std_s are the standard derivations
%   for closeness and similarity function.
%
%   Reference
%   ---------
%   This implemtation is based on the original paper 'Bilateral Filtering 
%   for Gray and Color Images' published by C. Tomasi and R. Manduchi 
%   (Proceedings of the 1998 IEEE International Conference on Computer Vision).
%
%   (c) Christopher Rohkohl
%       christopher@oneder.de
%       http://www.oneder.de
 
    % create gaussian closeness function
    hs    = (size-1)/2;
    [x y] = meshgrid(-hs(2):hs(2),-hs(1):hs(1));
    H_c   = exp(-(x.*x + y.*y)/(2*std_c*std_c));
 
    % perform filtering
    FA = nlfilter(A, size, @simfunc, std_s, H_c);
end
 
% adaptive similarity function
function V = simfunc(B, std, H_c);
    center = floor((size(B)+1)/2);
    sim = reshape(exp(-(B(:) - repmat(B(center(1), center(2)), numel(B), 1))...
                                    .^2 / (2*std*std)), size(B)) .* H_c;
    ssum = sum(sim(:));
    if (ssum ~= 0)
        sim = sim ./ sum(sim(:));
    end
    V = sum(sum(sim.*B));
end