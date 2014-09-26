This repository contains the analysis scripts and data associated with:
Chandler CH, Chari S, Tack D, Dworkin I. 2014. Causes and consequences of genetic background effects illuminated by integrative genomic analysis. Genetics. 196(4):1321-36. doi: 10.1534/genetics.113.159426 

This repository contains all analysis scripts and some raw/intermediate/final data files. All remaining original raw data are deposited at the [NCBI SRA](http://www.ncbi.nlm.nih.gov/Traces/sra/?study=SRP035237), and remaining intermediate files can be reconstructed from the raw data using the scripts included in this package. See readme for details.

To repeat analyses (see below for more description of directory structure and files):

I. Read mapping and genome re-sequencing

-- Download and decompress raw sequence reads for the SAM/ORE parental genomes in sequencing_scripts/reseq_data/sam_ore/, the backcross lines in sequencing_scripts/reseq_data/bc/, and the Drosophila reference genome in sequencing_scripts/reseq_data/reference/

-- Run sequencing_scripts/read_mapping/index_reference.sh to create bwa and samtools index files for the Drosophila reference genome

-- Run sequencing_scripts/read_mapping/map_*.sh to align reads to the reference genome; this will put intermediate data files in sequencing_scripts/temp_alignment_data/ and final alignments in sequencing_scripts/bam_alignments/

-- Run sequencing_scripts/snp_calling/call_all_snps.sh followed by sequencing_scripts/snp_calling/generate_pileup.sh, and subsequently generate_allele_freq_plot_data.sh. This will call SNPs, placing the results in sequencing_scripts/snp_data (as a VCF file and a pileup file), and will use custom scripts to calculate frequencies of the SAM and ORE alleles for the backcross line across the genome

-- Run sequencing_scripts/consensus_calling/call_*_consensus.sh to generate SAM and ORE genome sequences from the Drosophila reference genome and the SNP calls from earlier. Results will be placed in sequencing_scripts/consensus_sequences

II. Binding site analyses

-- Run sequencing_scripts/binding_sites_custom_meme/generate_expanded_resequened_cds.sh to generate fasta files with coding sequences + flanking regions for resequenced SAM and ORE genomes

-- Run sequencing_scripts/binding_sites_custom_meme/generate_motif_weights.sh to generate parameter values for binding site predictions

-- Run sequencing_scripts/binding_sites_custom_meme/motifscan_*_custom_meme.sh to run SD binding site predictions

-- Run sequencing_scripts/binding_sites_custom_meme/generate_motifscan_pvals.sh to run statistical analyses on binding site predictions for later

-- Run sequencing_scripts/targeted_sd_binding_site_analysis/make_pileup.sh to get pileup file for cut and sal SD binding sites

III. Gene expression

-- Run expression_scripts/microarray_analyses/do_microarray_analysis.sh to generate results for the Illumina and DGRC expression arrays

-- Run expression_scripts/dge/generate_ore_cds.sh, expression_scripts/dge/generate_sam_cds.sh, and finally expression_scripts/dge/generate_combined_cds.sh to create fasta files of the SAM and ORE coding sequences using the Drosophila reference genome and the SNP calls from earlier

-- Run expression_scripts/dge/index_cds.sh to index those files for read mapping

-- Run expression_scripts/dge/map_fc*.sh to map the DGE reads

-- Run expression_scripts/dge/pre_process_fc_mapping_results.sh, followed by expression_scripts/dge/combine_quantitative_data.sh and then expression_scripts/dge/remove_false_allele_specifics.sh to generate counts for downstream DGE analysis

-- Run expression_scripts/dge_r_analysis/do_dge_r_analysis.sh to complete DGE analyses (identify differentially expressed genes, etc.)

IV. Deletion mapping

-- Run src_scripts/do_src_analysis.sh

V. Putting everything together

-- Look for overlap/common signal across datasets: run integrating_datasets/allelic_imbalance/allelic_imbalance_to_arrays.R and integrating_datasets/allelic_imbalance/allelic_imbalance_to_binding.R; and integrating_datasets/compare_expression/overlap_tests.sh

-- Use the gene lists in integrating_datasets/go_analyses/ to perform GO analyses using external, web-based tools

-- Make plots of backcross allele frequencies: run integrating_datasets/plot_backcrosses/generate_and_plot_bc_bins.R

-- Make integrated plots summarizing across datasets: run integrating_datasets/plot_all/generate_bc_freqs.R followed by integrating_datasets/plot_all/compile_and_plot_with_direction.R



#####################################################################################
#####################################################################################

Data and scripts are organized into hierarchical folders structured as follows. In some cases, large files containing raw or processed data have been omitted from the actual Dryad data package, but can be regenerated using the included scripts and the raw data files that are available elsewhere (NCBI SRA).


1. sequencing_scripts/:

	-reseq_data/: (most data files are too large for Dryad but can be downloaded from NCBI SRA)
		-merge_samples.sh: combines raw data files for read mapping
		-bc/:
			-contains raw sequence reads for backcross lines
		-sam_ore/:
			-contains raw sequence reads for SAM & ORE
		-reference/:
			-contains reference genome

	-read_mapping/:
		-index_reference.sh: indexes reference genome sequences for bwa and samtools
		-map_*.sh: maps reads to reference genome

	-temp_alignment_data/:
		-contains temporary data files from read alignment (too big to be included in Dryad package but can be re-generated using scripts)

	-bam_alignments/:
		-contains BAM files from read mapping (too big to be included in Dryad package)

	-snp_calling/:
		-call_all_snps.sh: does actual SNP calling using samtools
		-generate_pileup.sh: generates a samtools pileup file at SNP sites
		-generate_allele_freq_plot_data.sh & .py: generates allele frequency data for backcross lines for binning and plotting
		
	-snp_data/:
		-contains raw SNP data (mostly too big for Dryad; Dryad data package does include VCF file)
	
	-consensus_calling/:
		-scripts that call consensus sequences and generate fasta files for the SAM and ORE genomes

	-consensus_sequences/:
		-contains consensus genome sequences for SAM & ORE (too big for Dryad)

	-targeted_sd_binding_analysis/:
		-contains scripts and data for known SD binding sites in cut and salm 

	-binding_sites_custom_meme/:
		-generate_expanded_resequenced_cds.py & .sh: generates transcript files for use in binding site prediction
		-sd.wtmx: position-weight matrix for SD binding site generated from meme
		-generate_motif_weights.sh: generates motif weights for use in MotifScan
		-motifscan_*_custom_meme.sh: does the actual motif scanning; generates LLR scores for SAM & ORE
		-calculate_base_comp.py, generate_motifscan_pvals.sh & .R: generates p-values from LLR scores generated earlier
		-SAM_ORE_binding_results_bygene_custom_meme.csv: the results of all of this
		
2. expression_scripts/:

	-microarray_data/:
		-contains data from previously published microarray experiments, along with processed summary results
	
	-microarray_analyses/:
		-contains scripts that generate the summary results of microarray experiments from the raw data

	-flybase/:
		-contains some data processed from flybase, and the script used to process the original flybase data

	-resequenced_cds/:
		-contains coding sequences for ORE & SAM genomes for mapping of DGE data; generated using the scripts above (too big for Dryad)

	-dge_data/:
		-contains raw and processed DGE reads and compiled summary results (SAM files are omitted from Dryad)
		 Note: FC3 and FC4 are sd[E3], while FC5 and FC6 are WT

	-dge/:
		-contains scripts to perform DGE read mapping and generate quantitative data (including scripts to generate coding sequences for SAM & ORE)
	
	-dge_temp_alignment/:
		-contains intermediate alignment files for DGE read mapping (not included in Dryad)
	
	-dge_r_analysis/:
		-contains scripts to analyze DGE experiments. Includes allelic_imbalance_rxn_norms.R, which generates PublicationAIReactionNorms.pdf, which is a figure in the paper
		

3. src_scripts/:

	-contains data and scripts for analyzing the deletion modifier screen (Chari & Dworkin 2013 PLoS Genetics)
	
4. integrating_datasets/:

	-functions/:
		-contains some R functions that are used in multiple analyses

	-allelic_imbalance/:
		-contains scripts that test for overlap between allelic imbalance significant genes, and binding predictions and genes significant in expression datasets

	-compare_expression/:
		-contains scripts that look for a congruent signal between expression datasets and binding site predictions


	-go_analyses/:
		-contains datasets of significant genes and reference gene lists for GO analyses

	-plot_backcrosses/:
		-generate_and_plot_mock_BC.R: generates the idealized introgression outcome plot
		-generate_and_plot_bc_bins.R: generates the plots of actual introgression

	-plot_all/:
		-generate_bc_freqs.R: generates allele frequencies in sliding windows for backcross lines for integrated plots (freq_data_10kbwindows_2kbsteps.csv)
		-compile_and_plot_with_direction.R: makes the integrated "everything" plots (contained in subdirectories)
		 also outputs a CSV file with all the compiled, integrated summary results (integrated_results_p.csv)
