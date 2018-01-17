#!/usr/bin/perl
use strict;
use warnings; 
use File::Copy;
use Cwd;

###################################################################
###################################################################
## 
## author: tanya
## 14.06.2012
##
## What it does?
## Takes a PSI-BLAST file as input 
## and calculates HSSP-values
## 
###################################################################
###################################################################


my $psiBlastDir=$ARGV[0]; # psi-blast file
#my $db=$ARGV[1]; #database name, i.e. extension of the psi-blast file (e.g. seq1_myDB.blastPsiOutTmp)


#my @psiBlastFiles = <$psiBlastDir/*_$db.blastPsiOutTmp>;
my @psiBlastFiles = <$psiBlastDir/*.blastPsiOutTmp>;

foreach my $psiBlastFile(@psiBlastFiles){
	convert2hssp($psiBlastFile);
#	print "printed output for $psiBlastFile\n";
}


###################################################################
###################################################################
# Subs
###################################################################
###################################################################
####################################################################
# calcs HSSP-Distance values from a regular BLAST or PSIBLAST file
####################################################################
sub convert2hssp {
my $fhin;
my $fhout;

my $file = shift;
# varieblaes for calculating HSSP-Distances for each alignment
my($line, $pi, $evalue, $num_identical, $num_gaps, $length, $hssp, $query, $hit, $query_ac, $hit_ac, $indxquery, $indxseq);

#check the file is a normal file and exists
die "file does not exist\n" if not -f $file;
open($fhin, '<', $file) || die ("could not open regular blast file $file in sub reg_blast_analysis\n");

#print "reading blast file $file for calculating hssp values\n";	

my $fileOut=$file . ".psiBlast2hssp";
print "$fileOut\n";

open ($fhout,">$fileOut") || die "could not open $fileOut\n";
print $fhout "###########################################################################\n";
print $fhout "# HIT - Name of a BLAST Hit\n";
print $fhout "# IDE - Percent identity in alignment\n";
print $fhout "# LALI - Length of alignment, excluding gaps\n";
print $fhout "# HSSP - HSSP-Distance using formula from B. Rost's HSSP-Paper 1999\n";
print $fhout "###########################################################################\n";
print $fhout "Hit\tIDE\tLALI\tEVAL\tHSSP\n";

my $isNewHit=0;
my $isConverged=0;

        while ($line = <$fhin>) {
reg_blast_1:
                if ($line =~ /^Query= (.+)/) {
                        $query = $1;
                        #my @values = split('#',$query);  # TODO change to correct separator
                        my @values = split(/\|/, $query);
                        $query_ac = $values[1];
                        $query = $query_ac;
			print $fhout "Query=$query\n";	
			$isConverged=0;
                }
		elsif ($line =~ /Results from round 3/ || $line =~ /CONVERGED!/ || $line =~ / No hits found /) {
   			$isConverged=1;
		}
#hit before whitespace                elsif ($line =~ /^[>](\S+) / && $isConverged==1) {
		elsif ($line =~ /^[>](.+)/ && $isConverged==1) {
                $hit = $1;
                my @hit_v = split('#',$hit);  # TODO change to correct separator
                #my @hit_v = split('|',$hit);
                $hit_ac = $hit_v[1];
                $hit = $hit_ac;
                $isNewHit=1;
	#	print "\thit=$hit\n";
reg_blast_3:
#                        next reg_blast_1 if ($query eq $hit);
reg_blast_4:
                        
#if($query ne $hit && $query_ac ne $hit_ac){  # this is the check if this is a selfhits
			if(1>0){  # this deactivates the check if this is a selfhits
                            #print "\n query -> ";
                            #print $query_ac;
                            #print "\n hit -> ";
                            #print $hit_ac;
                        while ($line = <$fhin>) {
                                #Score =  470 bits (1210), Expect = e-134
                                #Identities = 185/185 (100%), Positives = 185/185 (100%)
                                if (($line =~ /^\s*Score\s*=\s*[\d.]+\s*bits\s*\([\d.]+\),\s*Expect\s*=\s*([e\-.\d]+)/) && ($isNewHit==1)) {
                                        $line = <$fhin>;
					$evalue = $1;
                                        if ($line =~ /^\s*Identities\s*=\s*(\d+)\/(\d+)\s*\(\d+\%\),\s*Positives\s*=\s*(\d+)\/\d+\s*\(\d+\%\)(,\s*Gaps\s*=\s*(\d+)\/\d+\s*\(\d+\%\))?/) {
                                                $length = $2 - ($4 ? $5 : 0);
                                               # next reg_blast_4 if ($length < $options{'alilen'});
                                                $pi = ($1 / $length) * 100; # percent identities
                                                $hssp = sprintf("%.0f",hssp_dist($pi, $length));
                                                $pi = sprintf("%.1f",$pi);
                                               # $indxquery = $ids{$query}; $indxseq = $ids{$seq};
                                               # next line_tab if ($indxquery eq "" || $indxseq eq "");
                                               # $ali_data{$indxquery}{$indxseq} = $hssp if ((! $ali_data{$indxquery}{$indxseq}) || ($hssp > $ali_data{$indxquery}{$indxseq}));
                                        #       if ($hssp < $hsspthresh){
                                                        print $fhout "$hit\t$pi\t$length\t$evalue\t$hssp\n";
                                        #       }
                                        $isNewHit=0;
                                        }
                                        else {
                                             #   die( "ERROR: BLAST-format not correct:\nLine " . ($. - 1) . " is not followed by a line giving information about %identity\n");
						print( "ERROR: BLAST-format not correct:\nLine " . ($. - 1) . " is not followed by a line giving information about %identity\n");
                                        }
                                }
                              #  elsif ($line =~ /^[>](\S+)\|/) {
                              #   elsif ($line =~ /^[>](\S+) /) {
					elsif ($line =~ /^[>](.+)/) {
                                        $hit = $1;
                                        my @new_hit_v = split('#',$hit);  # TODO change to correct separator
                                        #my @new_hit_v = split(/\|/, $query);
                                        $hit_ac = $new_hit_v[1];
                                        $hit = $hit_ac;
                                        $isNewHit=1;
	#				print "\thit=$hit\n";
                                        goto reg_blast_3;
                                }
                                elsif ($line =~ /^Query= (.+)/) {
					$isConverged=0;
                                        goto reg_blast_1;
                                }
                        }
			}
                }
        }
if($isConverged ==0 ){
print "there was no round 3 word in $file\n";
print "program will exit!!!!\n\n";
exit 0;

}
close $fhout;
close $fhin;
return 1;
}

###################################################################
# calculates HSSP-Distance using the formula from Burkhard's
# HSSP-Paper 1999
# needs PIDE and length of alignment (LALI)
###################################################################
sub hssp_dist {

        my($pi) = shift; # PIDE = percent identity in alignment
        my($len) = shift; # LALI = length of alignment (minus gaps)

        if ($len <= 11) {
                return $pi - 100;
        }
        elsif ($len > 450) {
                return $pi - 19.5;
        }
        else {
                my($exp) = -0.32 * (1 + exp(- $len / 1000));
                return $pi - (480 * ($len ** $exp));
        }
}

####################################################################
# reads fa file and saves protein names (e.g. ABCD_HUMAN) in a hash
####################################################################
#sub getProteinsFromFa{

#my $fhin;
#	my $file = shift;
#	my %redprots; #hash of redundancy reduced proteins
	# check the file is a normal file and exists
#        die "file does not exist\n" if not -f $file;

#	open($fhin, '<', $file) || die ("could not open file $file\n");
#        print "reading blast file $file for calculating hssp values\n";
	
#	open ($fhout,">$fileOut") || die "could not open $fileOut\n";
#	while($line = <fhin>){
#		if($line = ~/^>/){
#			m/^>(.*?)\|/; #match shortest string between > and |
#    			$id=$1;
#			push (@{redprots{$id}});	
#		}
#	}
#	close $fhin;
#	return %redprots;
#}

