clear
clc

% ----------------------------------------------------------------------  %
% 
% Copyright (c) 2024 Aditi N Borkar, University of Nottingham
% All rights reserved.
% All codes, scripts and dataset in this package are distributed under the
% terms of the BSD 3-Clause License.
%
%
% This script removes 
%    - non-overlapping peaks between technical repeats and 
%    - overlapping peaks with gold reference.
% 
% 		This script takes as input 
% 		1. the gold reference peaks (gold.mat)
% 		2. the peak lists exported from SurfaceLab7 software containing 2
% 		columns of data as m/z values v/s intensity
% 		3. a simple nameslist file that lists the dataset names to be processed. Here,
% 			col1 = dataset name
% 			col2 = polarity
% 			col3 = total number of replicates
%	
% Expected Output: 
% 		data_peaks_1 struct that gives the processed m/z values for the RNA.
% 		Export it as a .mat file for further processing.
% 		e.g. output/data.mat 
% 			
% Notes:	
%		gold.mat contains a list of overlapping m/z values from two replicates collected
%		on blank gold sample substrate rinsed with PBS. 
%		
%		The path to input files required for this script may need to be modified depending
%		upon the user's directory structure and organisation of data
%		
%		The readtable and ismembertol functions used in this script are not yet implemented 
%		in GNU Octave, an open-source alternative to MATLAB. 
%		Octave users may need to consider using alternative functions 
%		or hard-code these functionalities in the current script for compatibility with Octave. 
%	
% ----------------------------------------------------------------------  %

% set tolerance
tol = 0.000001; % (1 ppm);

%import nameslist
nameslist=readtable('nameslist.txt', 'Format', '%s %s %s');

%import reference data
load('gold.mat');

for i=1:1:size(nameslist,1)
    
    disp(nameslist.Data_name{i});
    
    % determine polarity
    switch nameslist.polarity{i}
        case "neg"
            disp('Negative polarity data');
            ref_data = gold_neg_mz;
        case "pos"
            disp('Positive polarity data');
            ref_data = gold_pos_mz;
    end
    
    % determine number of replicates
    switch nameslist.replicates{i}
        case "1"
            disp('1 replicate processed')
            
            for j=1:1:1
                file=sprintf("orbisims_peaklist/%s-%s%d.txt", nameslist.Data_name{i},nameslist.polarity{i},j);
                spectrum = importdata(file);
                A{j} = spectrum(:,1);
            end
            
            data_peaks_1.([ nameslist.Data_name{i} ]) = A{1}(~ismembertol(A{1}, ref_data, tol)); % non-overlapping with reference.
            
        case "2"
            disp('2 replicates processed')
            
            for j=1:1:2
                file=sprintf("orbisims_peaklist/%s-%s%d.txt", nameslist.Data_name{i},nameslist.polarity{i},j);
                spectrum = importdata(file);
                A{j} = spectrum(:,1);
            end
            
            AB= A{1}(ismembertol(A{1} ,A{2} ,tol)); % find common peaks between repeats
            data_peaks_1.([ nameslist.Data_name{i} ]) = AB(~ismembertol(AB, ref_data, tol)); % non-overlapping with reference.
            
        case "3"
            disp('3 replicates processed')
            
            for j=1:1:2
                 file=sprintf("orbisims_peaklist/%s-%s%d.txt", nameslist.Data_name{i},nameslist.polarity{i},j);
                spectrum = importdata(file);
                A{j} = spectrum(:,1);
            end
            
            AB= A{1}(ismembertol(A{1} ,A{2} ,tol)); % find common peaks between repeats
            data_peaks_1.([ neg_names{i} ]) = AB(~ismembertol(AB, ref_data, tol)); % non-overlapping with reference.
            
            AC= A{1}(ismembertol(A{1} ,A{3} ,tol)); % find common peaks between repeats
            data_peaks_2.([ neg_names{i} ]) = AC(~ismembertol(AC, ref_data, tol)); % non-overlapping with reference.
            
            BC= A{2}(ismembertol(A{2} ,A{3} ,tol)); % find common peaks between repeats
            data_peaks_3.([ neg_names{i} ]) = BC(~ismembertol(BC, ref_data, tol)); % non-overlapping with reference.
            
        otherwise
            disp('Error: A maximum of 3 replicates can be processed')
    end
    
    clearvars file data A AB AC BC;
    
end

clearvars i j gold_neg_mz gold_pos_mz ref_data spectrum nameslist tol;