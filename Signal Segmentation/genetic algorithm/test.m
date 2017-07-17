total = [1 2 5 12 13 14 15 19 25 26 27 30]
actualSegments = [];
idx = 1;
flag = 0;
umm = [];
while (idx < size(total,2))
    currVal = total(idx);
    nextVal = total(idx+1);
    if((nextVal - currVal) < 2.5)
        umm = [umm currVal];
        flag = 1;
    else
        if(flag == 1)
            flag = 0;
            umm = [umm currVal];
            actualSegments = [actualSegments mean(umm)];
            umm = [];
        else
            actualSegments = [actualSegments currVal];
        end
    end  
    idx = idx + 1;
end
if(flag == 1)
    umm = [umm nextVal];
    actualSegments = [actualSegments mean(umm)];
else
    actualSegments = [actualSegments nextVal];
end
disp(actualSegments);

