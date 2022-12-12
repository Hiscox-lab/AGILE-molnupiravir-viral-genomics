# Characterisation of SARS-CoV-2 genomic variations in response to molnupiravir treatment in the AGILE Phase IIa clinical trial.

I’ah Donovan-Banfield<sup>1,2</sup> , Rebekah Penrice-Randal<sup>1</sup>, Hannah Goldswain<sup>1</sup>, Aleksandra M. Rzeszutek<sup>3</sup>, Jack Pilgrim<sup>3</sup>, Katie Bullock<sup>4</sup>, Geoffrey Saunders<sup>5</sup>, Josh Northey<sup>5</sup>, Xiaofeng Dong<sup>1</sup>, Yan Ryan<sup>1</sup>, Helen Reynolds<sup>6</sup>, Michelle Tetlow<sup>6</sup>,  Lauren E. Walker<sup>6,7</sup>, Richard FitzGerald<sup>6,7</sup>, Colin Hale<sup>7</sup>, Rebecca Lyon<sup>7</sup>, Christie Woods<sup>7</sup>, Shazaad Ahmed<sup>8</sup>, Dennis Hadjiyiannakis<sup>9</sup>, Jimstan Periselneris<sup>10</sup>, Emma Knox<sup>5</sup>, Calley Middleton<sup>5</sup>, Lara Lavelle-Langham<sup>4</sup>, Victoria Shaw<sup>4,11</sup>, William Greenhalf<sup>4</sup>, Thomas Edwards<sup>12</sup>,  David G. Lalloo<sup>13</sup>,  Christopher J. Edwards<sup>14,15</sup>, Alistair C. Darby<sup>1,16</sup>, Miles W. Carroll<sup>2,17</sup>, Gareth Griffiths<sup>5</sup>, Saye Khoo<sup>6</sup>, Julian A. Hiscox<sup>*1,2,18</sup>, Thomas Fletcher<sup>*2,12</sup>.

1 Department of Infection Biology and Microbiomes, Institute of Infection, Veterinary and Ecological Sciences, University of Liverpool, UK.<br>
2 NIHR Health Protection Research Unit in Emerging and Zoonotic Infections, Liverpool, UK.<br>
3 Department of Evolution, Ecology and Behaviour, Institute of Infection, Veterinary and Ecological Sciences, University of Liverpool, UK.<br>
4 GCPLab Facility, Institute of Systems, Molecular and Integrative Biology, University of Liverpool, UK.<br>
5 Southampton Clinical Trials Unit, University of Southampton, UK.<br>
6 Department of Molecular and Clinical Pharmacology, Institute of Translational Medicine, University of Liverpool, UK.<br>
7 NIHR Royal Liverpool and Broadgreen Clinical Research Facility, Liverpool, UK.<br>
8 NIHR Manchester Clinical Research Facility, Manchester University NHS Foundation Trust, Manchester, UK.<br>
9 NIHR Lancashire Clinical Research Facility, Lancashire Teaching Hospitals NHS Foundation Trust, Preston, UK.<br>
10 NIHR Kings Clinical Research Facility, King’s College Hospital NHS Foundation Trust, London, UK.<br>
11 The Clinical Directorate, University of Liverpool, UK.<br>
12 Centre for Drugs and Diagnostics, Liverpool School of Tropical Medicine, Liverpool, UK.<br>
13 Liverpool School of Tropical Medicine, Liverpool, UK.<br>
14 Human Development and Health School, University of Southampton, Southampton, UK.<br>
15 NIHR Southampton Clinical Research Facility, University Hospital Southampton NHS Foundation Trust, Southampton, UK.<br>
16 NIHR Health Protection Research Unit in Gastrointestinal Infections, Liverpool, UK.<br>
17 Wellcome Centre for Human Genetics, Nuffield Department of Medicine, University of Oxford, UK<br>
18 A*STAR Infectious Diseases Laboratories (A*STAR ID Labs), Agency for Science, Technology and Research (A*STAR), Singapore.<br>

### Article now published in Nature Communications (Open Access): https://www.nature.com/articles/s41467-022-34839-9

## 1. Sequencing of longitudinal nasopharyngeal swab samples using the Nimagen amplicon protocol (https://www.nimagen.com/gfx/protcol-NimaGen-covid-wgs_v19.pdf)
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

## Access to raw fastq data 
Deposited in the National Center for Biotechnology Information (NCBI) Short Read Archive (SRA)
Accessible here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA854613

## Contact us
Correspondence to: I'ah Donovan-Banfield (i.donovan-banfield@liverpool.ac.uk), Tom Fletcher (tom.fletcher@lstmed.ac.uk) or Julian Hiscox (julianh@liverpool.ac.uk)

