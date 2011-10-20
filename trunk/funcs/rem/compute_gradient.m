function outArr = compute_gradient(inArr)

len = length(inArr);
outArr = zeros(len,1);
%remember, inArr = i - (i+1)
for i = 2:(len-1)
    outArr(i) = (abs(inArr(i-1)) + abs(inArr(i)))/2;
end