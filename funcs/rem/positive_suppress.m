function outArr = positive_suppress(inArr, ratio)

len = length(inArr);
outArr = zeros(len,1);

for i = 1:len
    outArr(i) = suppress(inArr(i), ratio);
end