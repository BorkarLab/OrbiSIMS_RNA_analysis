function [fragRNA] = create_fragmentSeq(fileID, fragLen)
sequence = fastaread(fileID);
fragNum = length(sequence.Sequence) - (fragLen);
if ( fragNum == 0 )
    fragNum = 1;
end
for i = 1:1:fragNum
    fragRNA(i) = extractBetween(sequence.Sequence, i, i+fragLen-1);
end

