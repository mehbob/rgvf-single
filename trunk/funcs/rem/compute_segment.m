function segInfo = compute_segment(forwardDiff,threshold)
 
    forwardDiff(abs(forwardDiff) < threshold) = 0;
    
    lastPos = find(forwardDiff~=0,1,'last');
    startPos = 0;
    segInfo = [];
    num = 1;
    while(startPos < lastPos)
        startPos = find(forwardDiff ~= 0,1,'first');
        polarity = sign(forwardDiff(startPos));
        forwardDiff(1:startPos-1) = forwardDiff(startPos);
        endPos = find(sign(forwardDiff)~=polarity,1,'first')-1;
        forwardDiff(1:endPos)=0;
        segInfo{num}=[polarity,startPos,endPos];
        num = num + 1;
        
        startPos = endPos + 1;
    end
