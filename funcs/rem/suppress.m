function resVal = suppress(inVal,ratio)

%inVal is a scalar, so does resVal

if(nargin < 2)
    ratio = 0.1;
end

resVal = inVal;

if(inVal > 0)
    resVal = ratio * inVal;
end
