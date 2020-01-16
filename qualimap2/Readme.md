# Qualimap

Qualimap is a platform-independent application written in Java and R that provides both a Graphical User Interface (GUI) and a command-line interface to facilitate the quality control of alignment sequencing data.

Shortly, Qualimap:

	Examines sequencing alignment data according to the features of the mapped reads and their genomic properties
	Povides an overall view of the data that helps to to the detect biases in the sequencing and/or mapping of the data and eases decision-making for further analysis.


The main features offered by Qualimap are:

	fast analysis across the reference genome of mapping coverage and nucleotide distribution;
	easy-to-interpret summary of the main properties of the alignment data;
	analysis of the reads mapped inside/outside of the regions defined in an annotation reference;
	computation and analysis of read counts obtained from intersting of read alignments with genomic features;
	analysis of the adequacy of the sequencing depth in RNA-seq experiments;
	support for multi-sample comparison for alignment data and counts data;
	clustering of epigenomic profiles.


Usage: qualimap [tool] [options]

Available tools:

	bamqc            Evaluate NGS mapping to a reference genome

CLI usage example:
dx run Qualimap -iin=accepted_hits.bam -icmd="qualimap bamqc" -y
