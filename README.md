Running:
 - java -jar DGraph_v2.jar

Example usage:
 1) Make distance graph from short peptides or already-aligned sequences
   - Input: Fasta file
   - Choose "Run DGraph"
   - Choose "Generate PD scores" and check "PD score for aligned sequences"
   - Onscreen instructions say how to adjust the figure generation.
   - Press 's' to generate a PDF + PNG of the output

 2) Make distance graph from unaligned sequences (Sliding window approach)
   - Input: Fasta file
   - Choose "Run DGraph"
   - Choose "Generate PD scores." (Optionally) change the window size. 
   - Onscreen instructions say how to adjust the figure generation.
   - Press 's' to generate a PDF + PNG of the output

 3) Make distance graph from unaligned sequences (Alignment approach)
  - Input: Fasta file
  - Run clustalw, i.e. from https://www.genome.jp/tools-bin/clustalw
  - Download the .aln file. (Optional) Copy+paste full clustalw output to file.
  - Run jalview and convert the .aln to a fasta of sequences with dashes
    - Note: To do this, open the .aln in jalview, Ctrl+A then Ctrl+C.
    - Then, in a text editor Ctrl+V to paste the Fasta of the alignment
  - Now continue with the example 1, using this new fasta file.

How to generate input files:
 - Dgraph supplies a number of utilities for getting input into a usable
   format.
 - "Coloring" file format: Readme TODO

Running in a script:
 - Run
   java -Dhelp -jar DGraph_v2.jar
   to see a list of batch processing options.

Build dependencies:
 - Depends on core.jar, itext.jar, and pdf.jar from the Processing project.
   - Depends on older versions of Processing, namely Processing-1.5.1.
   - See https://processing.org/
 - Release jars have the appropriate versions of these dependencies packaged
   in.
