function outArr = compute_forward_difference(inArr)

len = length(inArr);
outArr = zeros(len,1);

for i = 2: len-1
    outArr(i) = inArr(i) - inArr(i+1);
end