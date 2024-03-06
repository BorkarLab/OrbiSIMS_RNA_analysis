# OrbiSIMS_RNA_analysis
Custom MATLAB scripts and demo dataset for processing and analysing cryo-OrbiSIMS data generated for RNA systems.

Copyright Statement

	Copyright (c) 2024 Aditi N Borkar, University of Nottingham
	
	All codes, scripts and dataset in this package are distributed under the terms of the BSD 3-Clause License:

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
	(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
	HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
	EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

System Requirements

	Operating System:

		macOS: Ventura 13.6.1 (tested) or later
		Windows: 10 or later (equivalent)
		Linux: Ubuntu 18.04 or later (recommended)

	Software:

		MATLAB: 2021a (tested). 
		Earlier or later versions may work, but functionality may differ and may require adjustments to the code.
		
		Bash shell: or equivalent terminal emulator and scripting language
		Most Linux distributions and macOS come with a bash shell pre-installed. 
		Windows users may need to install a bash shell emulator, such as Git Bash or Windows Subsystem for Linux (WSL).
		The bash scripting components may require adjustments depending on the specific Linux distribution you are using.
		
	Hardware:

		No specific hardware requirements beyond those necessary to run MATLAB and the chosen operating system.

Installation guide

	a. MATLAB:

	Download and install MATLAB 2021a (or a compatible version) from the 
	MathWorks website: https://www.mathworks.com/help/install/ug/download-without-installing.html or
	your institutional distribution. Follow the on-screen instructions during the installation process.
	
	b. Bash Shell (if needed):

		1. Windows:

			Option 1: Git Bash:
				Download and install Git for Windows from https://git-scm.com/download/win.
				During installation, choose the option to "Use Git from the Windows command line" and "Install Git Bash"
			Option 2: Windows Subsystem for Linux (WSL):
				Enable WSL by following the instructions here: https://learn.microsoft.com/en-us/windows/wsl/
				Install a Linux distribution of your choice (e.g., Ubuntu) within WSL.
		
		2. Linux and macOS:

			Bash shell is typically pre-installed on these systems. 
			You may verify its presence by opening a terminal and typing bash --version.
			
	c. Script Files:

	Download the MATLAB script files to a desired location on your system.
	(Optional) Organize the scripts: Consider creating a dedicated folder for your project to keep the scripts and any related files organised.
	
	Permissions (if needed):
	
		No specific permission changes are typically required for MATLAB scripts.
		
		On Linux and macOS, ensure any bash scripts have executable permissions. 
		You can do this by opening a terminal, navigating to the script directory, and running the following command for each script:
		
		Bash
		chmod +x script_name.sh
	
		**Use code with caution.
		Replace script_name.sh with the actual filename of your script.
		
Instructions on running the data (DEMO)

	**You may need to modify certain paths or adjust configurations within the scripts for them to function correctly on your specific system. 
	**Depending on the complexity of your data and processing power of your computer, 
	the time to run these scripts could range from a few minutes to a few days.
	
	1. Pre-processing the peaks list
	
	script: technical_process.m
		This script removes non-overlapping peaks between technical repeats and 
		removes the peaks overlapping with gold reference.
	
		This script takes as input 
		1. the gold reference peaks (gold.mat), 
		2. the peak lists exported from SurfaceLab7 software and 
		3. a simple nameslist file that lists the dataset names to be processed. Here,
			col1 = dataset name
			col2 = polarity
			col3 = replicate number
			
	datafiles: 
		OrbiSIMS peaklist files containing m/z values v/s intensity as two columns of data
		e.g. TAR_100nm-neg1.txt and TAR_100nm-neg2.txt

	Running instructions:
	Open MATLAB
	Navigate to the directory containing your MATLAB script and input datafiles 
	Run the script using the "Run" button the ribbon header or by typing the script name and pressing Enter.
	
	Expected Output: 
		data_peaks_1 struct that gives the processed m/z values for the RNA.
		Export it as a .mat file for further processing.
		e.g. output/data.mat 
		
		
	2. Peak Assignments:
	
		script: peaks_match_5_parallel_functions.m
	
			This takes as input 
			1. the RNA sequence in fasta format, 
			2. input dataset of processed peak list, and
			3. table of RNA ion fragment types.
			These files are placed in the inputs folder
			
			Dependent functions for this script: 
		
			create_fragmentDB.m			
			parallel_assignments.m
			create_fragmentSeq.m		
			matchIon_parallel.m
			
			All dependent m files need to be in the current folder. 
		
		Running instructions:
		Open MATLAB 
		Navigate to the directory containing your MATLAB script and input datafiles 
		Open the script file.
		Change the sequence input, dataset name, dataset mat and ionfragments in the file as needed.
		Run the script using the "Run" button the ribbon header or by typing the script name and pressing Enter.
		
		Expected Output: 
			list of peaksmatch files that enlist each match between theoretical database
			created for an RNA fragment of length f and the experimental spectrum.
			f will range from 3 to N, where N is the length of the RNA investigated. 
		
	3. Assignment analysis
	
		script: peaks_analysis.m
		
			This script takes as input
			1. all the peaksmatch files output from peaks_match_5_parallel_functions.m file
			2. the RNA sequence in fasta format
			** These files are placed in the current folder
		
		Running instructions:
		Open MATLAB 
		Navigate to the directory containing your MATLAB script and input datafiles 
		Open the script file.
		Change the sequence input and dataset name in the parent script file as needed.
		Run the script using the "Run" button the ribbon header or by typing the script name and pressing Enter.
		
		Expected Output: 
			1. Statement in the MATLAB console enlisting the minimim RMSE value assignments, 
			the corresponding fragment length and redundancy.
			e.g. For the example dataset of TAR RNA, RMSE = 3.5816 ppm for window fragment = 28 and redundancy = 3.6333
			
			** If the redundancy is > 8.3, it might be worthwhile to choose the assignment file with redundancy
			value immediately below 8.3. This might not be with the minimum with an RMSE. 
			
			2. Output file "assignment_frequency_{dataset}.txt" that contains the following data:
				col1 = residue number
				col2 = assignment frequency in the spectrum
				col3 = zscore of the assignment freqeuncy.	

	4. Use frequencies for 2D structure prediction

		-	predict secondary structures of all the sequences : RNAstructure 
			1. Do without any restraints. https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/Predict1/Predict1.html
				- just input the RNA sequence in the corresponding box.
				- right click and download the CT file. 
			
			2. Do with restraints. https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/Predict1/Predict1.html
				- input the RNA sequence in the corresponding box.
				- upload residue number v/s zscore assingment frequecies (col1 and col3 data only) file in the SHAPE constraints box. Keep all other parameters constant.
				- right click and download the CT file, rename approriately.
		
			3. Convert the CT file to dot bracket notation. https://rna.urmc.rochester.edu/RNAstructureWeb/Servers/ct2dot/ct2dot.html
				- upload the CT file and choose Structure Number 1.
				- Note the output.
			
				Example For TAR: 
				>ENERGY = -54.1  tar
				GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGGGAACCCACUGCUUAAGCCUCAAAAAAGCUUGCCUUGAGUGCUUCAAGUAGUGUGUGCCCGUCUGUUGUGUGACUCUGGAAAAAAAAAAA
				(((.(((((((((((.(((((...((((......)))))))))))))))))))))))(((((((((((((((((...........))))).)))).))))))))....((.(((........)))...))...........
			
				Example for TAR_restrained:
				>ENERGY = -69.7  tar_restrained
				GGUCUCUCUGGUUAGACCAGAUCUGAGCCUGGGAGCUCUCUGGCUAACUAGGGAACCCACUGCUUAAGCCUCAAAAAAGCUUGCCUUGAGUGCUUCAAGUAGUGUGUGCCCGUCUGUUGUGUGACUCUGGAAAAAAAAAAA
				(((.(((((((((((.(((((...((((......))))))))))))))))))))))).......((((.((......)))))).((.((((((..(((.(((((......)).)))))).)).)))).))...........
		
			4. Compare the 2D structures: http://beagle.bio.uniroma2.it/index.php
				- Enter the results from above.
			
				Example for TAR
				%Identity for TAR and TAR_restrained = 52.48%
				%Similarity for TAR and TAR_restrained = 68.09%
				p-value = 0.0001, thus significant differences between the two predictions.

	5. Convert 2D to 3D on trRosettaRNA or RNAComposer webservers
			- Input the above secondary structures and the fasta sequence
			- Download the predicted PDB
			- Compare / Visualise using TMAlign or UCSF Chimera software. 
			
	

	
