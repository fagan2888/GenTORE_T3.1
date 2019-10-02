function [output1,output2] = F3_pattern(B,A)
    SIZE = length(B)-length(A)+1;       
    match = zeros(1,SIZE);
    for i = 1:SIZE
        match(i) = max(all(B(i:i-1+length(A))==A));
    end
    output1 = find(match == 1);
    idx = find(diff(output1)==1);
    output2 = output1;
    output2(:,idx+1) = [];
end