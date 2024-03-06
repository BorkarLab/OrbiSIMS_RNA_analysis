function [ionDB5, ionDB3, fragSize, fragSizeB] = create_fragmentDB(seq)
% Script to generate theoretical ion fragmentation patterns for
% all possible nucleotide combinations.


bases = ["A", "G", "C", "U"];
phosphate = 97.9769;
hydroxyl = 18.0105;
proton	= 1.0078;
mzstd = [ 345.21, 361.21, 321.18, 322.17 ];
mzbase = [ 135, 151, 111, 112 ];

% load standard mz
mz5 = importdata('inputs/ionfragments/5mz');
mz3 = importdata('inputs/ionfragments/3mz');

% if subtracting 5' Phosphate to 5' OH
mass_diff = 79.9663;
mass_diff_2 = 0.0432;

w = strlength(seq);
fragSize = 0;
%full fragments
for x = 1:1:w
    nuc = extractBetween(seq, x, x);
    nuc_num = find(bases == nuc);
    fragSize = fragSize + mzstd(nuc_num);
end
fragSize = fragSize - (w-1).*hydroxyl + (w.*2 - 1).*proton - mass_diff + mass_diff_2;
fragSizeB = fragSize - (w-1).*hydroxyl + (w.*2 - 1).*proton - mzbase(nuc_num) - mass_diff;

% a - d ions
ionDB5 = NaN(w,size(mz5,2));
for x = 1:1:w-1
    nuc = extractBetween(seq, x, x);
    nuc_num = find(bases == nuc);
    
    for y = 1:1:size(mz5,2)
        
        if(x == 1)
            ionDB5(x,y) = mz5(nuc_num, y);
        else
            ionDB5(x,y) = mz5(nuc_num, y) + ionDB5(x-1,size(mz5,2)) - phosphate;
        end
    end
end
ionDB5 = ionDB5 - mass_diff;

% w - z ions
ionDB3 = NaN(w,size(mz3,2));
for x = w:-1:2
    nuc = extractBetween(seq, x, x);
    nuc_num = find(bases == nuc);
    
    for y = size(mz3,2):-1:1
        if(x == w)
            ionDB3(x,y) = mz3(nuc_num, y);
        else
            ionDB3(x,y) = mz3(nuc_num, y) + ionDB3(x+1,1) - hydroxyl;
        end
    end
end




