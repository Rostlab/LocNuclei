# LocNuclei

LocNuclei predicts 1) in which of the 13 subnuclear compartments a nuclear protein is located and 2) whether a nuclear protein is also located in another compartment outside the nucleus (traveller proteins).
LocNuclei combines homology-based inference if possbile with a de nove prediction using a string-based Profile Kernel SVM otherwise.

## Requirements
LocNuclei is a Python project. To run LocNuclei, Python 3 has to be installed locally. To run LocNuclei for a protein, the protein sequence as fasta format and the BLAST profile for this protein is needed.

## How to use it
LocNuclei can easily be called from the command line. The call needs to include at least a folder to fasta-files, a folder to profile-files and the path to an output file in which the results will be written.
One can use different parameters to modify the programs behavior.

**Positional arguments**
**fasta_folder**: Folder with protein sequences in fasta format. Every file may only contain one sequence.
**blast_folder**: Folder with BLAST profiles for the fasta files. BLAST files need to have the same name in front of the suffix as the fasta files (e.g. Q9XLZ3.fasta <-> Q9XLZ3.profile)
**output_file**: Path to a file in which the results shall be written

**Optional arguments**
**-h, --help**: Shows a help message and explains the parameters
**--fasta_suffix**: Suffix of files in given fasta-folder (default: "\*.fasta")
**--blast_suffix**: Suffix of files in given BLAST folder (default:"\*.profile")
**--temp_folder**: Folder to work in, will be automatically deleted afterwards. If not given, a temporary directory will b e created.
**-t, --traveller**: Predict nuclear travelling proteins instead of sub-nuclear localization
**-b, --only_blast**: Run only homology based inference
**-d, --debug**: Toggles clean up of temporary files off, i.e. no files will be deleted that were created during prediction
**-v, --verbose**: Toggles verbose mode on

## Output
The output of LocNuclei contains one line per protein and per line 4 columns. The 1st column is the protein Id, followed by all predicted location classes in the format "LocA. LocB.". The 3rd column contains the source (b = blast, s = svm) and the last columns contains the reliability index (RI). For homology based inference, there is one RI derived from the percentage pairwise sequence identity of the BLAST hit. For de novo prediction, there is one RI for every predicted subnuclear compartment.

## Example
The folder **example/** contains three fasta and profiles files to try out LocNuclei. To run the example, a simple call inthe folder where **locnuclei.py** is located is sufficient:

`python locnuclei.py example/ example/ result.out -t` 



