# AGILE CST-2 Viral Genomics project
 - Phase IIa trial investigating molnupiravir (EIDD-2801) efficacy

## 1 Sequencing longitudinal nasopharyngeal swab samples using the Nimagen amplicon protocol (https://www.nimagen.com/gfx/protcol-NimaGen-covid-wgs_v19.pdf)
    1. 536 samples in total (180 patients days 1, 3 & 5, except 3 patients who missed >= 1 visits) 
    2. Amplicon library prep and QC carried out by Jack Pilgrim and Alex Rzeszutek 
    3. Illumina sequencing conducted by Centre for Genomics Research
    4. Demultiplexing and initial data processing conducted by Richard Gregory and Sam Haldenby (CGR Informatics team)
    5. Raw fastq files processed with the EasySeq pipeline (v0.9; https://github.com/JordyCoolen/easyseq_covid19) using the default parameters 
    - mutation frequency: >=0.5
    - QUAL:	>=20
    - Minimum sequence depth: 10
    - Maximum %N: 0.09999
    6. Pangolin lineage assignments, minor variant analysis and all data visualisation done by I'ah Donovan-Banfield
        - 177 patients have day 1 lineage data.


Sample details in sample_info directory: 
1. AGILE-nimagen-full-180_final-metadata.csv (includes Panoglin lineage info, N count and treatment allocation information)
2. SRA accession (PRJNA854613) metadata files 
    - SRA-metadata_AGILE-library-prep.xlsx
    - SRA_submission_metadata_AGILE_sctu.xlsx



## Scripts available `/home/idb75/scripts`:
1. diversitools_loop_opt.sh - XD/RPR/HG (updated to include Syn_NonSyn_parse_aa_V3_loop_IDB.pl by IDB)
2. diversiutils_aa_details.pl - XD
3. Syn_NonSyn_parse_aa_V3_loop_IDB.pl - IDB
4. nucleotide_analysis.R
5. amino-acid_analysis.R



Contributor initials: IDB - I'ah Donovan-Banfield; RPR - Rebee Penrice-Randal; HG - Hannah Goldswain; XD - Xiaofeng Dong


# Workflow for Paper 
To add later 

