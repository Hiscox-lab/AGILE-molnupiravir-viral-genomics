# AGILE CST-2 Viral Genomics project

## 1. Sequencing longitudinal nasopharyngeal swab samples using the Nimagen amplicon protocol (https://www.nimagen.com/gfx/protcol-NimaGen-covid-wgs_v19.pdf)
    1. 536 samples in total (180 patients days 1, 3 & 5, except 3 patients who missed >= 1 visits) 
    2. Amplicon library prep and QC carried out by Jack Pilgrim and Alex Rzeszutek 
    3. Illumina sequencing conducted by Centre for Genomics Research (CGR)
    4. Demultiplexing and initial data processing conducted by Richard Gregory and Sam Haldenby (CGR Informatics team)
    5. Raw fastq files processed with the EasySeq pipeline (v0.9; https://github.com/JordyCoolen/easyseq_covid19) using the default parameters 
    - mutation frequency: >=0.5
    - QUAL:	>=20
    - Minimum sequence depth: 10
    6. Pangolin lineage assignments, minor variant analysis and all data visualisation done by I'ah Donovan-Banfield
        - 177 patients have day 1 lineage data.

Sample details in sample_info directory: 
1. AGILE-nimagen-full-180_final-metadata.csv (for stratifying data in R analysis)
2. SRA accession (PRJNA854613) metadata files 
    - SRA-metadata_AGILE-library-prep.xlsx
    - SRA_submission_metadata_AGILE_sctu.xlsx

## Workflow 
### 1. Run raw fastq files through the EasySeq pipeline 
   - v0.9 used, accessed at https://github.com/JordyCoolen/easyseq_covid19
### 2. Concatenate consensus.fasta files 
   Run: 
   - pangolin 
   - nextclade
   - facount (calculates N count)
### 3. Check coverage breadth
   - each patient must have all 3 samples with genome breadth of 90% at 200X 
   - move all passed bam files to new directory
### 4. Generate metadata table with sample / trial information 
   Includes:
   - Kit id (unique sample identifier)
   - Treatment group 
   - Visit day 
   - N count 
   - Lineage information 
   
### 5. Run DiversiTools on all passed samples
   - requires the primer trimmed bam file (called sample-name_L001.final.bam in EasySeq outputs)
   - Run the diversitools_loop_opt.sh script (uses diversiutils_aa_details.pl and Syn_NonSyn_parse_aa_V3_loop_IDB.pl)
### 6. Nucleotide analysis
   - Conducted in R with the entropy.txt DiversiTools output files
### 7. Amino acid analysis 
   - Conducted in R with the AA.txt DiversiTools output files

### Scripts available:
1. diversitools_loop_opt.sh - XD/RPR/HG (updated to include Syn_NonSyn_parse_aa_V3_loop_IDB.pl by IDB)
2. diversiutils_aa_details.pl - XD
3. Syn_NonSyn_parse_aa_V3_loop_IDB.pl - IDB
4. nucleotide_analysis.R
5. amino-acid_analysis.R

Contributor initials: IDB - I'ah Donovan-Banfield; RPR - Rebee Penrice-Randal; HG - Hannah Goldswain; XD - Xiaofeng Dong



