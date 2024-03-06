% ----------------------------------------------------------------------  %
% 
% Copyright (c) 2024 Aditi N Borkar, University of Nottingham
% All rights reserved.
% All codes, scripts and dataset in this package are distributed under the
% terms of the BSD 3-Clause License.
%
%
% This script creates a theoretical dataset of all ions arising 
% from a given RNA sequence and matches it with the experimental peaklist.
% This script takes as input
%       1. the RNA sequence in fasta format,
%       2. input dataset of processed peak list, and
%       3. table of RNA ion fragment types.
%    ** These files are placed in the inputs folder
% 
% Dependent functions for this script:
%       1. create_fragmentDB.m
%       2. parallel_assignments.m
%       3. create_fragmentSeq.m
%       4. matchIon_parallel.m
%    ** All dependent .m files need to be in the current folder.
% 
% Change the sequence input, dataset name, dataset mat and ionfragments 
% in the parent script file as needed.
% 
% Expected Output: 
% 		A number of peaksmatch files that enlist each match between theoretical database
% 		created for an RNA fragment of length f and the experimental spectrum.
% 		f will range from 3 to N, where N is the length of the RNA investigated. 
%
% Notes:	
%		The path to input files required for this script may need to be modified depending
%		upon the user's directory structure and organisation of data
%		
%		fastaread function used in this script would require BioInformatics toolbox in MATLAB
%		and the bioinfo 0.2.0 package in GNU Octave.
%		
%		The writematrix and ismembertol functions used in this script and its dependent functions
%		are not yet implemented in GNU Octave, an open-source alternative to MATLAB. 
%		Octave users may need to consider using alternative functions 
%		or hard-code these functionalities in the current script for compatibility with Octave. 
%		
% 		
% ----------------------------------------------------------------------  %

clear
clc

%constants
tol = 1E06; % tolerance 1 ppm
proton = 1.0078;
iontypes5 = [ "a-B", "a", "b", "c", "d"];
iontypes3 = ['w', 'x', 'y', 'z'];

ionloss = ["neutral", "OH", "H2O", "HPO3", "H3PO4" ];
ionlossVal = [0, 17.0027, 18.0105, 79.9663, 97.9769 ];


% RNA sequence file
seqFile = "inputs/tar.fa";

% import peak list
% already processed to remove overlaps with reference and
% common pekas between repeat scans
dataset = {"TAR_100nm"};
load('inputs/data.mat', 'data_peaks_1');
spectrum = data_peaks_1.( [ dataset{1} ]);
spectrum_tol = spectrum / tol * 2;

% determine window ranges
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

for w = 1:1:size(wndw,1)
    fragRNA = create_fragmentSeq(seqFile, wndw(w));
    
    fileID = sprintf('peaksmatch_%s.%s.txt',dataset{1}, num2str(wndw(w)));
    header = ["Window", "Seq", "theory_m/z", "theory_m/z", "Peak", "ion", "deletions", "Score" ];
    writematrix(header,fileID,'Delimiter','tab');
    
    % suppress E notation
    format longG
    
    parfor i = 1:1:size(fragRNA,2) % for each window of RNA sequence
        
        [m1, m2, fragSize, fragSizeB] = create_fragmentDB(fragRNA(i));             % copy standard ion fragments to local variables.
        
        tic
        % assignment of intact segments
       
        x = 1;
        y = 1;
        chk = "segment";
        parallel_assignments(wndw(w), spectrum, spectrum_tol, i, fragRNA(i), fragSize, fileID, chk, x, y);
         
        % assignment of ions
        chk = "5ion";
        for j = 1:1:size(m1,1)-1 % for each nucleotide in fragment
            for k = 1:1:size(m1,2) % for each iontype
                parallel_assignments(wndw(w), spectrum, spectrum_tol, i, fragRNA(i), m1(j,k), fileID, chk, k, j);
            end
        end
        
        chk = "3ion"
        n = size(m2,1);
        for j = 2:1:size(m2,1)
            n = n - 1;
            for k = 1:1:size(m2,2)
                parallel_assignments(wndw(w), spectrum, spectrum_tol, i, fragRNA(i), m2(j,k), fileID, chk, k, n);
            end
        end
        
        elapsed_time(w,i) = toc;
    end
end