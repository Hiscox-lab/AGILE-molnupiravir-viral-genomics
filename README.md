# Characterisation of SARS-CoV-2 genomic variations in response to molnupiravir treatment in the AGILE Phase IIa clinical trial.

I’ah Donovan-Banfield1,2, Rebekah Penrice-Randal1, Hannah Goldswain1, Aleksandra M. Rzeszutek3, Jack Pilgrim3, Katie Bullock4, Geoffrey Saunders5, Josh Northey5, Xiaofeng Dong1, Yan Ryan1, Helen Reynolds6, Michelle Tetlow6,  Lauren E. Walker6,7, Richard FitzGerald6,7, Colin Hale7, Rebecca Lyon7, Christie Woods7, Shazaad Ahmed8, Dennis Hadjiyiannakis9, Jimstan Periselneris10, Emma Knox5, Calley Middleton5, Lara Lavelle-Langham4, Victoria Shaw4,11, William Greenhalf4, Thomas Edwards12,  David G. Lalloo13,  Christopher J. Edwards14,15, Alistair C. Darby1,16, Miles W. Carroll2,17, Gareth Griffiths5, Saye Khoo6, Julian A. Hiscox*1,2,18, Thomas Fletcher*2,12.

1 Department of Infection Biology and Microbiomes, Institute of Infection, Veterinary and Ecological Sciences, University of Liverpool, UK.
2 NIHR Health Protection Research Unit in Emerging and Zoonotic Infections, Liverpool, UK.
3 Department of Evolution, Ecology and Behaviour, Institute of Infection, Veterinary and Ecological Sciences, University of Liverpool, UK.
4 GCPLab Facility, Institute of Systems, Molecular and Integrative Biology, University of Liverpool, UK.
5 Southampton Clinical Trials Unit, University of Southampton, UK
6 Department of Molecular and Clinical Pharmacology, Institute of Translational Medicine, University of Liverpool, UK
7 NIHR Royal Liverpool and Broadgreen Clinical Research Facility, Liverpool, UK
8 NIHR Manchester Clinical Research Facility, Manchester University NHS Foundation Trust, Manchester, UK
9 NIHR Lancashire Clinical Research Facility, Lancashire Teaching Hospitals NHS Foundation Trust, Preston, UK
10 NIHR Kings Clinical Research Facility, King’s College Hospital NHS Foundation Trust, London, UK
11 The Clinical Directorate, University of Liverpool, UK 
12 Centre for Drugs and Diagnostics, Liverpool School of Tropical Medicine, Liverpool, UK
13 Liverpool School of Tropical Medicine, Liverpool, UK
14 Human Development and Health School, University of Southampton, Southampton, UK
15 NIHR Southampton Clinical Research Facility, University Hospital Southampton NHS Foundation Trust, Southampton, UK
16 NIHR Health Protection Research Unit in Gastrointestinal Infections, Liverpool, UK.
17 Wellcome Centre for Human Genetics, Nuffield Department of Medicine, University of Oxford, UK
18 A*STAR Infectious Diseases Laboratories (A*STAR ID Labs), Agency for Science, Technology and Research (A*STAR), Singapore.



raw fastq data available from the National Center for Biotechnology Information (NCBI) Short Read Archive (SRA) (Project Accession Number PRJNA854613)

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
4. nucleotide_analysis.R - IDB
5. amino-acid_analysis.R - IDB

Contributor initials: IDB - I'ah Donovan-Banfield; RPR - Rebee Penrice-Randal; HG - Hannah Goldswain; XD - Xiaofeng Dong

Correspondence to: I'ah Donovan-Banfield (i.donovan-banfield@liverpool.ac.uk), Tom Fletcher (tom.fletcher@lstmed.ac.uk) or Julian Hiscox (julianh@liverpool.ac.uk)

