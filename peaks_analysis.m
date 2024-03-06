% ----------------------------------------------------------------------  %
% 
% Copyright (c) 2024 Aditi N Borkar, University of Nottingham
% All rights reserved.
% All codes, scripts and dataset in this package are distributed under the
% terms of the BSD 3-Clause License.
%
%
% This script finds the assignment with minimum RMSE and/or redundancy
% and converts the m/z assignments to peak assignment frequency for each residue in the RNA seqeunce. 
% 
% This script takes as input
%       1. all the peaksmatch files output from peaks_match_5_parallel_functions.m
%       file
%       2. the RNA sequence in fasta format
%    ** These files are placed in the current folder
% 
% Change the sequence input and dataset name in the parent script file as needed.
% 
% Expected Output: 
% 			1. Statement in the MATLAB console enlisting the minimim RMSE value assignments, 
% 			the corresponding fragment length and redundancy.
% 			e.g. For the example dataset of TAR RNA, RMSE = 3.5816 ppm for window fragment = 28 and redundancy = 3.6333
% 			
% 			** If the redundancy is > 8.3, it might be worthwhile to choose the assignment file with redundancy
% 			value immediately below 8.3. This might not be with the minimum with an RMSE. 
% 			
% 			2. Output file "assignment_frequency_{dataset}.txt" that contains the following data:
% 				col1 = residue number
% 				col2 = assignment frequency in the spectrum
% 				col3 = zscore of the assignment freqeuncy.	
%
% Notes:	
%		The path to input files required for this script may need to be modified depending
%		upon the user's directory structure and organisation of data
%		
%		fastaread function used in this script would require BioInformatics toolbox in MATLAB
%		and the bioinfo 0.2.0 package in GNU Octave.
%		
%		The readtable, writematrix and ismembertol functions used in this script 
%		are not yet implemented in GNU Octave, an open-source alternative to MATLAB. 
%		Octave users may need to consider using alternative functions 
%		or hard-code these functionalities in the current script for compatibility with Octave. 
%		
% 
% ----------------------------------------------------------------------  %

clear
clc

% load all the files

% RNA sequence file
seqFile = "tar.fa";
dataset = {"TAR_100nm"};

sequence = fastaread(seqFile);
seqLen = length(sequence.Sequence);
clearvars sequence;

p = 1;
wndw = zeros(seqLen,1);
while ( ( floor(seqLen / p ) >= 3 ) )
    wndw(p) = floor(seqLen / p );
    p = p + 1;
end
wndw(p) = floor ( (floor(seqLen/2) + seqLen)/2 );
p = p + 1;
wndw(p) = seqLen - 1;
p = p + 1;
wndw(p:end) = [];
temp = unique(wndw);
wndw = temp;
clearvars temp;

% find assignment with minimum RMSE and/or redundancy
for i=1:1:size(wndw,1)
    file=sprintf("peaksmatch_%s.%d.txt", dataset{1}, wndw(i));
    assignments = readtable(file);
    
     % save window size
     data.window{i} = wndw(i);
     
     % initiate statistics
     data.spectrum_redundancy{i} = 0;
     data.rmse{i} = 0;
     data.rmse_2{i} = 0;
     
     if ( size(assignments,1) ~= 0 )
         
         %calculate redundancy
         [GC,GR] = groupcounts(assignments.Var5);
         data.spectrum_redundancy{i} = sum(GC)/size(GC,1);
         clearvars GC, GR;
         
         %calculate RMSE
         % counts
         [GC,GR] = groupcounts(assignments.Var1);
         % error
         error = zeros(size(GR,1));
         for j=1:1:size(GR,1)
             peaks_theory = assignments.Var4(assignments.Var1(:) == GR(j));
             peaks_exp = assignments.Var5(assignments.Var1(:) == GR(j));
             error(j) = sqrt ( sum((peaks_theory - peaks_exp).*(peaks_theory - peaks_exp)) / GC(j));
         end
         data.rmse{i} = sum(error(:))/size(GR,1);
         data.rmse_2{i} =  data.rmse{i}/size(GR,1)*10^06;
         data.assigned{i} = size(GR,1);
         clearvars GC GR error;
     end
end

A1 = cell2mat(data.rmse_2);
rmse_min = min(A1(A1 > 0 ));

A2 = cell2mat(data.window);
rmse_wndw = A2(A1 == rmse_min);

A3 = cell2mat(data.spectrum_redundancy);
rmse_redundancy = A3(A1 == rmse_min);

str = sprintf("RMSE = %8.4f ppm for window = %d and redundancy = %8.4f\n"...
    , rmse_min, rmse_wndw, rmse_redundancy);
disp(str);

clearvars assignments A1 A2 A3 error GC GR i j p ;

% convert assignments to frequency
file=sprintf("peaksmatch_%s.%d.txt", dataset{1},rmse_wndw);
assignments = readtable(file);

frequency = zeros(seqLen,1);
ion5 = [ "a-B", "a", "b", "c", "d" ];
for i=1:1:size(assignments,1)
    resid = 0;
    iontype = [];
    ion=regexp(assignments.Var6{i},'_','split');
    if ( size(ion,2) == 1 )
        iontype = 'M';
        resid = assignments.Var1(i) + rmse_wndw - 1;
    else
        iontype = extract(string(ion{1}),1);
        resnum = str2num(extractAfter(string(ion{1}),1));
        if ( contains(iontype, ion5) )
            resid = assignments.Var1(i) + resnum - 1;
        else
            resid = assignments.Var1(i) + rmse_wndw - resnum;
        end
    end
    frequency(resid) = frequency(resid) + 1;
end
frequency_norm = zscore(frequency);

 fileID = sprintf('assignment_frequency_%s.%s.txt',dataset{1}, num2str(rmse_wndw));
 header = ["ResNum", "Frequency_count", "Frequency_zscore"];
 writematrix(header,fileID,'Delimiter','tab');
 
 for i=1:1:seqLen
     str = sprintf("%8d%10d%12.4f",i,frequency(i), frequency_norm(i));
     writematrix(str,fileID,'WriteMode','append','Delimiter','tab');
 end
   
