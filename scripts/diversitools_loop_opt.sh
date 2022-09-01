#!/bin/bash
#make directory with all your bams and bai files in

while getopts  "a:b:r:c:n:l:s:e:h" c 
do 
    case $c in 
	a) ANALYSIS=$OPTARG ;;
	r) REFERENCE=$OPTARG ;;
	c) CODING=$OPTARG ;;
	b) BAMS=$OPTARG ;;
	n) NONSYN=$OPTARG ;;
	l) LOADING_LIST=$OPTARG ;;
	s) SEDCOM=$OPTARG ;;
	e) EXT=$OPTARG ;;
	h) echo " usage sh diversitools_loop_opts.sh "
echo " 		-a [output analysis directory ie 20210728_output] " 
echo "          -b [absolute path to directory with all bam and bai files for DT analysis]"
echo "		-r [absolute path to reference fasta file] "
echo "		-c [absolute path to coding regions txt file] "
echo "		-n [absolute path to syn_nonsynscript] "
echo "		-l [file name of loading list including extension -contains file names of all input files, for nonsyn script]"
echo "		-s [sed command of common string at end of all your file names including the separating underscore i.e 's/_L002.final.bam//' "
echo "		-e [extension of file names that you remove in sed command, without the sed bits so just <_L002.final.bam> ***has to be identical to what you remove in SEDCOM***"
exit 0 ;;
    esac 
done

#make analysis directory with directory within that has the sym link to all the bam and bai files you'll need
echo creating output directory
mkdir -p ${ANALYSIS}
touch ${ANALYSIS}/DTlog_$(date -u).txt
echo "script ran on $(date -u)" > ${ANALYSIS}/DTlog_$(date -u).txt
echo "The command used was sh diversitools_loop_opts.sh -a $ANALYSIS -r $REFERENCE -c $CODING -b $BAMS -n $NONSYN -l $LOADING_LIST -s $SEDCOM -e $EXT" >> ${ANALYSIS}/DTlog_$(date -u).txt
echo entering bam directory
cd ${BAMS}

for i in $(ls *bam | sed ${SEDCOM} | uniq)
do
        echo "processing sample $i" >> ${ANALYSIS}/DTlog_$(date -u).txt
		/home/xiaofeng/sources/perl_5.32.1/bin/perl ~/projects/agile_mpv_seq_project/scripts/diversiutils_aa_details.pl -bam ${i}${EXT} -ref $REFERENCE -orfs $CODING -stub ${ANALYSIS}/${i}

done

#make loading list for syn_nonsyn script
echo making loading list
ls *bam | sed ${SEDCOM} > ${LOADING_LIST}

#cd ../
echo running syn_nonsyn script
cd ${ANALYSIS}
#run non syn script
# use script that has optional arguments from ~/projects/hiscox-artic-pipeline/scripts/Syn_NonSyn_parse_aa_V3_loop.pl 
# usage path/to/script/ ARG0 <loading_list_selected.tsv> ARG1 < path to location of AA.txt and AA_parse.txt > 
/home/xiaofeng/sources/perl_5.32.1/bin/perl $NONSYN ${BAMS}/${LOADING_LIST} $ANALYSIS
echo "diversitools and syn_nonsyn scripts completed on $(date -u)" >> ${ANALYSIS}/DTlog_$(date -u).txt
