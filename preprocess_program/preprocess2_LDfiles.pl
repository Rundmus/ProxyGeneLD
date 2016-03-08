use strict;
use warnings;

# ------------------------------------------------------------------------------
# Function: trim out low r_squared LD and compact files
# ------------------------------------------------------------------------------
#
# Usage   : >perl preprocess2_LDfiles.pl <Arguments>
#
# ------------------------------------------------------------------------------

my $startTime = time;
my $chkptTime = $startTime;

use File::Basename;
use File::Spec;

my $tolerance = 10000;	# Gene closed to SNP was defined as gene within 10kbp

sub showUsage {
	print <<END;
Usage: >perl preprocess2_LDfiles.pl (snpF) (gF) (pth) (pS) [th=.8]

     snpF - hapmapSNP file from UCSC table browser
     gF   - the file that was created by preprocess1_GeneIDnSym.pl
     pth  - the path includes HapMap LD data files
     pS   - population symbol used in HapMap LD data files (e.g. CEU)
     th   - threshold of r2 to trim out lower values than (default 0.8)
END
	exit;
}
showUsage if ( scalar(@ARGV) < 4 );

my $snpFN  = shift(@ARGV);
my $geneFN = shift(@ARGV);

my $path = shift(@ARGV);
my $popSym = shift(@ARGV);
my $threshold = shift(@ARGV) || 0.8;
my $outFNprefix = 'hapmap' . $popSym . '_GeneID_r2_' . $threshold . '_chr';

print "Threshold of r_squared to trim out is " . $threshold . "\.\n";

###
### parse information about reference genes
###
my %posBin = ();
{
	print "parsing " . fileparse($geneFN) . "\n";
	
	open(GENE, $geneFN) or die "Can't open $geneFN: $!";
	<GENE>; 			# ignore the title line
	while(<GENE>) {
		chomp;
		my ($gid, undef, @d) = split /\t/;		# split the attributes
		
		#		 @d     0     	1   	2    	3     	4      	5    	    6       7    	 	8      		9        	10
		#GeneID	Symbol	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	id	name2	cdsStartStat	cdsEndStat	exonFrames
		#79501	OR4F5	NM_0..	chr1	+		58953	59871	58953		59871	1			58953,		59871,		0	OR4F5	cmpl	cmpl	0,
		
		my $millionStart = int(($d[3] - $tolerance) / 1000000);
		my $millionEnd   = int(($d[4] + $tolerance) / 1000000);
		
		my @exSts = split (/\,/, $d[8]);
		my @exEnds   = split (/\,/, $d[9]);
		
		my %indDatum = (
			'gid'     => $gid,
			'strand'  => $d[2],
			'txSt'    => $d[3]  , 'txEnd'   => $d[4],
			'cdsSt'   => $d[5]  , 'cdsEnd'  => $d[6],
			'exCount' => $d[7],
			'exSts'   => \@exSts, 'exEnds'  => \@exEnds
		);
		
		$d[1] =~ s/chr//;
		for (my $i = $millionStart; $i <= $millionEnd; $i++) {
			my $bin = $d[1] . '_' . $i;			# bin : chr.num._millionNum.
			                           			#     e.g. chr10 54213135 -> bin is 10_54
			####################################
			push (@{$posBin{$bin}}, \%indDatum);
			####################################
		}
		
		if (abs($d[4] - $d[3]) > 2350000) {
			print "$d[0] is longer than 2.35 million bps.\n";
		}
	}
	showElapsedTimeSince($chkptTime);
	$chkptTime = time;
}

###
### parse location info of HapMap SNPs
###
my %hapmapSNP = ();				# 'chr number' -> array of SNPs
{
	print "parsing " . fileparse($snpFN) . "\n";
	open(SNP, $snpFN) or die "Can't open $snpFN: $!";
	<SNP>;			# ignore title line
	
	while(<SNP>) {
		chomp;
		
		#bin	chrom	chromStart	chromEnd	name	score	strand	observed	allele1	homoCount1	allele2	homoCount2	heteroCount
	
		my (undef, $chr, undef, $pos, $rsID) = split /\t/;
		$rsID = lc($rsID);
		
		$chr =~ s/chr//;		# remove chr prefix to leave only chr number
		my @snpSet = ($rsID, $pos);
		
		#####################################
		push (@{$hapmapSNP{$chr}}, \@snpSet);
		#####################################
	}
	close SNP;
}

###
### parse LD info files
###

{
	$| = 1;			# forces a flush after every print to show progress
	use List::Util qw(max min);
	
	foreach my $chr ( 1..22, 'X' ) {
		## file names
		my $inFN = File::Spec->catfile($path, 'ld_chr' . $chr . '_' . $popSym . '.txt');
		my $outFN = File::Spec->catfile($path, $outFNprefix . $chr);
		
		showElapsedTimeSince($chkptTime);
		$chkptTime = time;
		print "Working on Chr$chr...";
		
		open(IN, $inFN) or die "Can't open $inFN: $!";
		
		my %snpInLD = ();
		while(<IN>) {
			#554636 558390 CEU rs9629043 rs11497407 1.0 0.0010 0.04 5
			
			my (undef, undef, undef, $rs1, $rs2, undef, $rsqr) = split;
			
			next if($rsqr < $threshold);
			
			$rs1 = lc($rs1);
			$rs2 = lc($rs2);
			
			if(defined $snpInLD{$rs1}) {
				$snpInLD{$rs1} .= "\t" . $rs2 . ' ' . $rsqr;
			}else {
				$snpInLD{$rs1} = $rs2 . ' ' . $rsqr;
			}
		}
		close IN;
		
		open(OUT, ">$outFN") or die "Can't open $outFN: $!";
		print OUT join("\t", 'rs#', 'chr', 'pos', 'GeneIDs', 'from_Genes', 'from_CDSs', "from_closest_exon_5\'(3\'ss)", "from_5\'ss", "SNPs_in_LD r2") . "\n";
		
		foreach my $snp (@{$hapmapSNP{$chr}}) {
			my ($rsID, $pos) = @$snp;
			
			my $bin = $chr . '_' . int($pos / 1000000);
			print OUT join("\t", $rsID, $chr, $pos) . "\t";
		
			if ( exists $posBin{$bin} ) {
				my %g2txs = ();
				foreach my $tx ( @{$posBin{$bin}} ) {		# $tx : hash of transcript all info
					# within $tolerance bps from gene
					if ((($tx->{'txSt'}-$tolerance) < $pos) && ($pos <= ($tx->{'txEnd'}+$tolerance))) {
						push @{$g2txs{$tx->{'gid'}}}, $tx;
					}
				}
				if ( keys %g2txs ) {
					my @sortedGeneIds = sort keys %g2txs;
					print OUT join("\,", @sortedGeneIds);
					
					my (@dGs, @dCdss, @dEx5, @dEx3) = ();
					foreach my $gid (@sortedGeneIds) {
						# for each gene
						my ($dGMin, $dCdsMin, $dEx5Min, $dEx3Min) = (1000000) x 4; # shortest distances among all transcripts of the gene
						my $upstream = 0;			# if SNP is located in upstream of CDS
						foreach my $tx (@{$g2txs{$gid}}) {		# for each transcript of the gene
							my $distTxPt = $tx->{'txSt'} + 1 - $pos;
							my $distTxQt = $pos - $tx->{'txEnd'};
							if ($distTxPt <= 0 && $distTxQt <= 0) {
								$dGMin = 0;			# within transcript, dGMin (distance from gene) = 0
							}else {
								$dGMin = min(abs($distTxPt), abs($distTxQt), $dGMin);
							}
							
							my $dCdsPt = $tx->{'cdsSt'} + 1 - $pos;
							$upstream = 1 if ($dCdsPt > 0 && $tx->{'strand'} eq "+");
							my $dCdsQt = $pos - $tx->{'cdsEnd'};
							$upstream = 1 if ($dCdsQt > 0 && $tx->{'strand'} eq "-");
							if ($dCdsPt <= 0 && $dCdsQt <= 0) {
								$dCdsMin = 0;		# within CDS, $dCdsMin = 0
							}else {
								$dCdsMin = min(abs($dCdsPt), abs($dCdsQt), $dCdsMin);
							}
							
							############################
							## searching closest exon ##
							############################
							next if $dGMin != 0;
							
							if($tx->{'strand'} eq "-") {
								for (my $i = 0;$i < $tx->{'exCount'};$i++) {
									if ( ($pos - $tx->{'exEnds'}->[$i]) <= 0 ) {		# if (distance to exon end <= 0)
										my $distExSt = $tx->{'exSts'}->[$i] + 1 - $pos;
										if ( $distExSt <= 0 ) {							# if (distance to exon start <= 0)
											$dEx3Min = max($distExSt, -1* abs($dEx3Min));		# closest to splice site
											$dEx5Min = max($pos - $tx->{'exEnds'}->[$i], -1* abs($dEx5Min));
											last;
										}else {
											$dEx3Min = min($distExSt, $dEx3Min);
											$dEx5Min = min($pos - $tx->{'exEnds'}->[$i-1], $dEx5Min);
											last;
										}
									}
								}
							}else {
								for (my $i = 0;$i < $tx->{'exCount'};$i++) {
									if ( ($pos - $tx->{'exEnds'}->[$i]) <= 0 ) {		# if (distance to exon end <= 0)
										my $distExSt = $tx->{'exSts'}->[$i] + 1 - $pos;
										if ( $distExSt <= 0 ) {							# if (distance to exon start <= 0)
											$dEx5Min = max($distExSt, -1* abs($dEx5Min));		# closest to splice site
											$dEx3Min = max($pos - $tx->{'exEnds'}->[$i], -1* abs($dEx3Min));
											last;
										}else {
											$dEx5Min = min($distExSt, $dEx5Min);
											$dEx3Min = min($pos - $tx->{'exEnds'}->[$i-1], $dEx3Min);
											last;
										}
									}
								}
							}
						}
						
						# still for each gene
						if($upstream && $dGMin != 0) { $dGMin *= -1; }		# upstream -> negative
						push @dGs, $dGMin;
						if($upstream && $dCdsMin != 0) { $dCdsMin *= -1; }	# upstream -> negative
						push @dCdss, $dCdsMin;
						
						if($dGMin != 0) { $dEx5Min = $dEx3Min = ""; }
						
						push @dEx5, $dEx5Min;
						push @dEx3, $dEx3Min;
					}
					print OUT "\t" . join("\t", join("\,", @dGs), join("\,", @dCdss), join("\,", @dEx5), join("\,", @dEx3));
					
				}else {
					print OUT "i", "\t" x 4;
				}
			}else {
				print OUT "i", "\t" x 4;
			}
			print OUT "\t";
			print OUT $snpInLD{$rsID} if(exists $snpInLD{$rsID});
			print OUT "\n";
		}
		close OUT;
	}
	showElapsedTimeSince($chkptTime);
	$| = 0;		# return forces a flush after every print
}

showElapsedTimeSince($startTime);

sub showElapsedTimeSince {	# show elapsed time
	my $since = shift;
	my $prefix = shift;
	$prefix = 'time elapsed : ' unless defined $prefix;
	my $suffix = shift;
	$suffix = "\n" unless defined $suffix;

	if (!defined $since) {
		print "There seems to be an error in calling subroutine \'showElapsedTime\'.";
		return;
	}
	
	my $elapsedSec = time - $since;
	my $elapsedMin = int(($elapsedSec + 0.1)/ 60);
	$elapsedSec -= $elapsedMin*60;
	my $elapsedHour = int(($elapsedMin + 0.1)/ 60);
	$elapsedMin -= $elapsedHour*60;
	
	print $prefix;
	if ( $elapsedHour ) {
		printf "%d h %d m %d s" , $elapsedHour, $elapsedMin , $elapsedSec;
	}elsif ( $elapsedMin) {
		printf "%d m %d s" , $elapsedMin , $elapsedSec;
	}else {
		printf "%d sec." , time - $since;
	}
	print $suffix;
}