function outArr = stack_refine(inArr, segInfo)

len = length(inArr);

segSize = size(segInfo);
segNum = segSize(2);

positiveArea = 0;

for i = 1 :  segNum
    if segInfo{i}(1) == 1 
        positiveArea = positiveArea + ...
                        sum(inArr( segInfo{i}(2) : segInfo{i}(3) ));
    else if positiveArea > 0
            j = segInfo{i}(2);
            while(j <= segInfo{i}(3) && positiveArea > 0)
                tmp = inArr(j); % remember tmp is negative
                assert(tmp < 0, 'tmp is assumed to be negative!');
                inArr(j) = min(tmp + positiveArea, 0);
                positiveArea = max(positiveArea + tmp, 0);
                j = j + 1;
            end
        end
    end
end

outArr = inArr;


