#!/usr/bin/perl -w
# usage Syn_NonSyn_parse_aa_V3_loop_IDB.pl ARG0 ARG1; ARG0= loading_list_selected.tsv; ARG1= path to location of AA.txt and AA_parse.txt

open (INDEX, "$ARGV[0]");
@indexs=<INDEX>;
chomp @indexs;
close INDEX;


foreach $index(@indexs) {
%hash=();
open(DATA,"$ARGV[1]/$index\_aa_details.txt");
while (<DATA>) {
    chomp;
    @each=split(/\t/);
    push (@{$hash{"$each[0]\_$each[4]"}}, "$each[5]\_$each[6]\_$each[7]\_$each[8]\_$each[1]");
}
close DATA;

open (NEWDETAIL,">$index\_AA_parse.txt");
print NEWDETAIL "Sample\tChr\tProtein\tAAPosition\tRefAA\tRefSite\tRefCodon\tCntNonSyn\tCntSyn\tNbStop\tAAcoverage\n";
    
open(AA, "$ARGV[1]/$index\_AA.txt");
while (<AA>) {
    unless (/^Sample\tChr/) {
        @eachindex=split(/\t/);
        $Chr=$eachindex[1];
        $protein=$eachindex[2];
        $AAPosition=$eachindex[3];
        $RefAA=$eachindex[4];
        $RefSite=$eachindex[5];
        $RefCodon=$eachindex[6];
       
        if ($eachindex[4]=~/\<NA\>/) {
            print NEWDETAIL "$index\t$Chr\t$protein\t$AAPosition\t0\t0\t0\t0\t0\t0\t0\n";
        }else{
            $CntNonSyn=0;
            $CntSyn=0;
            $stop=0;
            $AAcoverage=0;
            
            foreach $compare(@{$hash{"$AAPosition\_$protein"}}) {
                #print "$AAPosition\_$eachindex[0]\t$compare\n";
                $AAcoverage++;
            
                @eachcompare=split(/\_/,$compare);

                if ($eachcompare[1] eq $eachcompare[3] && $eachcompare[0] ne $eachcompare[2] && $eachcompare[2] ne "NNN") {
                    $CntSyn++;
                }
                     
                if ($eachcompare[1] ne $eachcompare[3]) {
                    $CntNonSyn++;
                }
                  
                if ($eachcompare[1]!~/\*/ && $eachcompare[3] =~/\*/ && $eachcompare[2] ne "NNN") {
                    $stop++;
                }
            }
            print NEWDETAIL "$index\t$Chr\t$protein\t$AAPosition\t$RefAA\t$RefSite\t$RefCodon\t$CntNonSyn\t$CntSyn\t$stop\t$AAcoverage\n";
        }
    }
}
close AA;
close NEWDETAIL;
}
