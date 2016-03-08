#!/usr/bin/perl
use strict;
use warnings;

# ------------------------------------------------------------------------------
# Function: Reading through SNP list, it searches Entrez GeneID of the gene
# 		enough close to each SNP or SNP in LD with it. 
#		It creates new file that contains adjusted gene-wide allelic P-value.
# ------------------------------------------------------------------------------
#
# Usage   : >perl proxyGeneLD <Arguments>
#
# INPUT files
#	1: \PATH_proxyGeneLD
#	2: (query snp list file)
#
# OUTPUT files
#	1: filename : "query snp list file name" . "_GeneLD" . threshold of r2
#
#		"\t", geneId  symbol  position  rsId  unadjP  adjP  SNP#nearGene  SNP#incl.LD  adjSNP#  length
#
# ------------------------------------------------------------------------------

my $startTime = time;
my $chkptTime = $startTime;

use File::Basename;

sub showUsage {
	print <<END;
Usage: >perl proxyGeneLD (qry) (rs) (P) (d) (pS) [th] [r2] [5'] [3']

     qry - file of query snp list
     rs  - column number containing rs num of SNP
     P   - column number of P-value
     d   - delimiter (space, tab, comma)
     pS  - population symbol used in HapMap LD data files (e.g. CEU)
     th  - threshold to trim out in preprocess2_LDfiles.pl (default = 0.8)
     r2  - threshold of r2 (def = 0.8)
     5'  - 5' promoter region length of gene (def = 1000)
     3'  - possible 3'UTR region length (def = 0)
END
	exit;
}
showUsage if ( scalar(@ARGV) < 5 );

## Read PATH file to find where the LD files are

my $pathFile = 'PATH_proxyGeneLD';
open(PTHS, $pathFile) or die "Can't open $pathFile: $!";
my @paths = ();
while(<PTHS>) {
	chomp;
	next if($_ =~ /^\#/);
	push @paths, $_;
}

## file names
my $geneFN = $paths[0];			# for pre-mRNA length calculation

my $ldPath  = $paths[1];

## Arguments
my $qrySnpFN 	 = shift(@ARGV);
my $colRs 		 = shift(@ARGV) - 1;
my $colPvalue 	 = shift(@ARGV) - 1;
my $paramDelimit = lc(shift(@ARGV));
my $popSym       = uc(shift(@ARGV));
my $thLDfiles    = shift(@ARGV) || 0.8;
my $r2threshold  = shift(@ARGV) || 0.8;
my $promoLen	 = abs(shift(@ARGV) || 1000) * (-1);
my $utr3Len		 = abs(shift(@ARGV) || 0);

if($thLDfiles > $r2threshold) {
	print qq(\nThe threshold used in trimming by preprocess2 must be ) .
		qq(smaller than the threshold for defining proxies, ) .
		qq(that is \(th\)<=\(r2\)\.\n);
	exit;
}

my $ldFNprefix = 'hapmap'. $popSym .'_GeneID_r2_'. $thLDfiles .'_chr';
		# where the files of SNP and LD info are

my $outFN  = $qrySnpFN . "_LD" . $r2threshold;		# output file name
my $nonHapMapFN = $qrySnpFN;

my $delimiter = '';
use Switch;
switch ($paramDelimit) {
	case "tab" { $delimiter = '\t' }
	case "space" { $delimiter = ' ' }
	case "comma" { $delimiter = ',' }
	else { showUsage }
}

###############################
## parse query snp list file ##
###############################

print "Reading " . basename($qrySnpFN) . "\n";

################################################################################

my %inQrySnp = %{createHash_fromFile($qrySnpFN, $delimiter, $colRs, $colPvalue)};

## %inQrySnp : Was the SNP in the Query Snp list
##           = ( rs# => P-value, )
################################################################################

showElapsedTimeSince($chkptTime);
$chkptTime = time;

my $countSnp = scalar(keys %inQrySnp);

#####################################################################################
### parse information about reference genes in order to calculate pre-mRNA length ###
#####################################################################################

my %gSPL = ();			# geneID -> symbol, position, average length of all transcripts
{
	print "Parsing " . basename($geneFN) . "\n";
	open(GENE, $geneFN) or die "Can't open $geneFN: $!";
	<GENE>; 			# ignore the title line
	
	my %gLensTx = ();
	while(<GENE>) {
		chomp;
		my ($gid, $symbol, undef, $chr, undef, $txSt, $txEnd, undef) = split (/\t/, $_, 8);
		
		#0     	1   	2    	3     	4      	5    	6       7    	 	8      	9        	10
		#GeneID	Symbol	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	id	name2	cdsStartStat	cdsEndStat	exonFrames
		#79501	OR4F5	NM_004	chr1	+		58953	59871	58953		59871	1	 		58953,		59871,		0	OR4F5	cmpl	cmpl	0,
		
		next if($chr =~ /\_/);		# exclude e.g. chr5_random, chr22_h2_hap1
		
		push @{$gLensTx{$gid}}, ($txEnd - $txSt);
		if (exists $gSPL{$gid}) {
			my (undef, $txStEnd) = split(/:/, $gSPL{$gid}->[1]);
			my ($txStBefore, $txEndBefore) = split(/-/, $txStEnd);
			$txSt = $txStBefore if ($txSt > $txStBefore);
			$txEnd = $txEndBefore if ($txEnd > $txEndBefore);
		}
		@{$gSPL{$gid}} = ($symbol, $chr . ':' . $txSt . '-' . $txEnd);
	}
	close GENE;

	use List::Util qw(sum);
	
	foreach my $gid (keys %gLensTx) {
		$gSPL{$gid}->[2] = sum(@{$gLensTx{$gid}}) / @{$gLensTx{$gid}};
	}
	
	showElapsedTimeSince($chkptTime);
	$chkptTime = time;
}

$| = 1;			# forces a flush after every print to show progress

my %prnOut = ();
my %snpCountOf = ();

use File::Spec;

## for each chromosome
#foreach my $i (4) {#debug
foreach my $i ( 1..22, 'X' ) {
	
	###
	### parse hapmapCEU_GeneID_LD_chr** files
	### store GeneID and location info of HapMap SNPs
	###

	## file name that contains info of HapMap SNPs
	my $snpFN = File::Spec->catfile($ldPath, $ldFNprefix . $i);
	
	my @outSnps = ();
	my %genesOf = ();
	
	## show progress
	print "Working on Chr$i...";

	################################################################################
		
	my ($rLDUp, $rSnpList) = parse_SNPfile_readLD($snpFN);

	## %{$rLDUp} : Link between SNPs in LD
	##			 = ( rs# of SNP in downstream => rs# of SNP in upstream => undef, )
	## @{rSnpList} = (	\@(SNP rs num, geneIDs, distance from Genes), ) 
	################################################################################
	
	foreach my $rThisSnpData (@$rSnpList) {
		my ($thisSNP, $gid, $dGs) = @$rThisSnpData;
		###
		###
		###
		if (exists $inQrySnp{$thisSNP}) {		# Is this SNP one of the query SNP ?
			push @outSnps, $thisSNP;
			
			if ($gid eq 'i') {			# Does this SNP located in intergenic region?
				if (exists $rLDUp->{$thisSNP}) {		# Is there any SNP uptream in LD with this SNP?
					my %countIncrYes = ();
					my %countIncrNo = ();
					
					foreach my $eaLDUp ( keys %{$rLDUp->{$thisSNP}} ) {		# for each SNP in LD with this SNP
					
						## If any SNP in LD uptream was in any Gene, 
						## this SNP is considered to be in the gene.
						if (exists $genesOf{$eaLDUp}) {
							foreach my $eaGeneUp (keys %{$genesOf{$eaLDUp}}) {
								if (exists $inQrySnp{$eaLDUp}) {
									# if snp in LD upstream in gene are one of query snps
									# then ignore this snp which in intergenic region
									$countIncrNo{$eaGeneUp} = undef;
								}else {
									$genesOf{$thisSNP}{$eaGeneUp} = undef;
									$countIncrYes{$eaGeneUp} = undef;
								}
							}
						}
					}
					foreach my $key ( keys %countIncrNo ) {
						delete $countIncrYes{$key};
					}
					foreach my $key ( keys %countIncrYes ) {
						$snpCountOf{$key}{'count'} ++;
					}
				}
			
			## this query SNP located near to Gene(s)
			}else {
				
				my @gids = split /\,/, $gid;
				my @dGss = split /\,/, $dGs;
				
				my %thisSelfGids = ();
				for my $i (0 .. $#gids) {
					if (isClose2Gene($dGss[$i])) {
						my $thisGid = $gids[$i];
						$genesOf{$thisSNP}{$thisGid} = undef;
						$snpCountOf{$thisGid}{'count'} ++;
						$snpCountOf{$thisGid}{'pre'} = $snpCountOf{$thisGid}{'count'} - 1;
						$snpCountOf{$thisGid}{'loc'} ++;
						$thisSelfGids{$thisGid} = undef;
					}
				}
				if (exists $rLDUp->{$thisSNP}) {		# Is there any SNP uptream in LD with this SNP
					my %countIncrYes = ();
					my %countIncrNo = ();
					
					foreach my $eaLDUp ( keys %{$rLDUp->{$thisSNP}} ) {		# for each SNP in LD with this SNP
						if (exists $genesOf{$eaLDUp}) {
							foreach my $eaGeneUp (keys %{$genesOf{$eaLDUp}}) {
								if (exists $inQrySnp{$eaLDUp}) {
									# if snp in LD upstream in gene are one of query snps
									# then ignore this snp which in other genic region
									$countIncrNo{$eaGeneUp} = undef;
									
									## if a query SNP in LD with this was in gene, return back to previous value
									if (exists $thisSelfGids{$eaGeneUp}) {
										$snpCountOf{$eaGeneUp}{'count'} = $snpCountOf{$eaGeneUp}{'pre'};
									}
								}else {
									next if (exists $thisSelfGids{$eaGeneUp});
									$genesOf{$thisSNP}{$eaGeneUp} = undef;
									$countIncrYes{$eaGeneUp} = undef;
								}
							}
						}
						# if snp in LD upstream located out of genic region are one of query snps
						# then do NOT consider the snp belongs to the gene
						# eventhough this SNP in LD are in gene.
					}
					foreach my $key ( keys %countIncrNo ) {
						delete $countIncrYes{$key};
					}
					foreach my $key ( keys %countIncrYes ) {
						$snpCountOf{$key}{'count'} ++;
					}
				}
			}
			
		## this SNP is NOT in query list
		}else {
			
			if ($gid eq 'i') {		# intergenic region
				delete $rLDUp->{$thisSNP};
			}else {						# close to gene
				my @gids = split /\,/, $gid;
				my @dGss = split /\,/, $dGs;
				
				my %thisSelfGids = ();
				for my $i (0 .. $#gids) {
					if (isClose2Gene($dGss[$i])) {
						my $thisGid = $gids[$i];
						$genesOf{$thisSNP}{$thisGid} = undef;
						$thisSelfGids{$thisGid} = undef;
					}
				}
				
				if (keys(%thisSelfGids) && exists $rLDUp->{$thisSNP}) {		# Is there any SNP uptream in LD with this SNP
					foreach my $selfGid (keys %thisSelfGids) {
						my $hasNewQrySnp = 0;
#						my $test = '';#debug
						foreach my $eaLDUp ( keys %{$rLDUp->{$thisSNP}} ) {		# for each SNP in LD with this SNP
							if (exists $inQrySnp{$eaLDUp}) {
								unless (exists $genesOf{$eaLDUp}{$selfGid}) {
									$genesOf{$eaLDUp}{$selfGid} = undef;
									$hasNewQrySnp = 1;
#									$test = $eaLDUp;#debug
								}
							}
						}
						if($hasNewQrySnp) { 
#							print join("\t", $selfGid, $snpCountOf{$selfGid}{'count'},$thisSNP, $test) . "\n";#debug
							$snpCountOf{$selfGid}{'count'} ++; 
						}
					}
				}
			}
		}
	}
	
	my %genesInChr = ();
	
	foreach my $rsNum (@outSnps) {
		my $snpPvalue = delete $inQrySnp{$rsNum};
		
		next unless(exists $genesOf{$rsNum});
		
		foreach my $gid (keys %{$genesOf{$rsNum}}) {
			my $genePvalue = $snpPvalue * $snpCountOf{$gid}{'count'};
			
			if(! exists($prnOut{$gid}) || $genePvalue < $prnOut{$gid}->[2]) {
				@{$prnOut{$gid}}[0..2] = ($rsNum, $snpPvalue, $genePvalue);
			}
			$prnOut{$gid}->[3] ++;
			$genesInChr{$gid} = undef;
		}
	}
	
	foreach my $gidChr1 (keys %genesInChr) {
		foreach my $gidChr2 (keys %genesInChr) {
			next if($gidChr1 == $gidChr2);
			if($prnOut{$gidChr1}->[0] eq $prnOut{$gidChr2}->[0] ||
				exists($rLDUp->{$prnOut{$gidChr1}->[0]}->{$prnOut{$gidChr2}->[0]})) {
				$prnOut{$gidChr1}->[4]->{$gidChr2} = undef;
				$prnOut{$gidChr2}->[4]->{$gidChr1} = undef;
			}
		}
	}
				
	## show progress
	print "finished in ";
	showElapsedTimeSince($chkptTime, '');
	$chkptTime = time;
}

print "\n";
$| = 0;		# return forces a flush after every print

{
	$outFN .= '_' . scalar( keys %prnOut ) . '.out';
	my $outSigFN = $outFN . 'Trim';

	open(OUT,">$outFN") or die "Can't open $outFN: $!";
	open(UTSIG,">$outSigFN") or die "Can't open $outSigFN: $!";
	
	my %pushSig = ();
	my $countSig = 0;
	foreach my $gid (sort { $prnOut{$a}->[1] <=> $prnOut{$b}->[1] || $prnOut{$a}->[2] <=> $prnOut{$b}->[2] } keys %prnOut) {
		next if(exists $pushSig{$gid});
		foreach my $gidSig (sort { $prnOut{$a}->[2] <=> $prnOut{$b}->[2] || $prnOut{$a}->[1] <=> $prnOut{$b}->[1] } keys %{$prnOut{$gid}->[4]}) {
			next if($prnOut{$gidSig}->[0] eq $prnOut{$gid}->[0] || exists($pushSig{$gidSig}));
			$pushSig{$gidSig} = $countSig;
			$countSig ++;
		}
	}
	
	## %prnOut =         'geneId' =>                     (0       1         2                       3                                   4)
	print OUT join("\t", 'geneId', 'symbol', 'position', 'rsId', 'unadjP', 'adjP', 'SNP#nearGene', 'SNP#incl.LD', 'adjSNP#', 'length', 'geneInLD') . "\n";
	print UTSIG join("\t", 'geneId', 'symbol', 'position', 'rsId', 'unadjP', 'adjP', 'SNP#nearGene', 'SNP#incl.LD', 'adjSNP#', 'length', 'geneInLD') . "\n";
	
	foreach my $gid (sort { $prnOut{$a}->[2] <=> $prnOut{$b}->[2] || $prnOut{$a}->[1] <=> $prnOut{$b}->[1] } keys %prnOut) {
		print OUT join("\t", $gid, @{$gSPL{$gid}}[0, 1], @{$prnOut{$gid}}[0..2], 
			($snpCountOf{$gid}{'loc'} || 0), $prnOut{$gid}->[3], 
			$snpCountOf{$gid}{'count'}, $gSPL{$gid}->[2]). "\t";
		my @prnOutSig = ();
		foreach my $gidSig (sort { $prnOut{$a}->[2] <=> $prnOut{$b}->[2] || $prnOut{$a}->[1] <=> $prnOut{$b}->[1] } keys %{$prnOut{$gid}->[4]}) {
			push @prnOutSig,$gSPL{$gidSig}[0];
		}
		print OUT join(",", @prnOutSig) . "\n";
		
		next if(exists $pushSig{$gid});
		print UTSIG join("\t", $gid, @{$gSPL{$gid}}[0, 1], @{$prnOut{$gid}}[0..2], 
			($snpCountOf{$gid}{'loc'} || 0), $prnOut{$gid}->[3], 
			$snpCountOf{$gid}{'count'}, $gSPL{$gid}->[2], join(",", @prnOutSig)) . "\n";
	}
	print UTSIG "\n";
	foreach my $gid (sort { $prnOut{$a}->[2] <=> $prnOut{$b}->[2] || $prnOut{$a}->[1] <=> $prnOut{$b}->[1] } keys %pushSig) {
		my @prnOutSig = ();
		foreach my $gidSig (sort { $prnOut{$a}->[2] <=> $prnOut{$b}->[2] || $prnOut{$a}->[1] <=> $prnOut{$b}->[1] } keys %{$prnOut{$gid}->[4]}) {
			push @prnOutSig,$gSPL{$gidSig}[0];
		}
		print UTSIG join("\t", $gid, @{$gSPL{$gid}}[0, 1], @{$prnOut{$gid}}[0..2], 
			($snpCountOf{$gid}{'loc'} || 0), $prnOut{$gid}->[3], 
			$snpCountOf{$gid}{'count'}, $gSPL{$gid}->[2],
			join(",", @prnOutSig)) . "\n";
	}
	
	close UTSIG;
	close OUT;
}

{
	$nonHapMapFN .= '_' . scalar(keys %inQrySnp) . '.nonHapMap';
	open(NONHM,">$nonHapMapFN") or die "Can't open $nonHapMapFN: $!";
	
	foreach my $snp (sort {$inQrySnp{$a} <=> $inQrySnp{$b}} keys %inQrySnp) {
		print NONHM join("\t", $snp, $inQrySnp{$snp}) . "\n";
	}
	close NONHM;
}
print "There are ". scalar(keys %inQrySnp) ." non-HapMap SNPs in total $countSnp.\n";

showElapsedTimeSince($startTime);

#End of main

#####
# global variable used : $utr3Len, $promoLen

sub isClose2Gene {
	my ($distance2Gene) = @_;
	
	if ( $distance2Gene <= $utr3Len && $distance2Gene >= $promoLen ) {	# within transcript region of gene
		return 1;
	}else {
		return 0;
	}
}


################################
### Create Hash reading file ###
################################
# arguments : ( file name, delimiter, column of key, column of value)
# return : (refHash)
#	%{refHash} = ( lc key => value, )
#

sub createHash_fromFile {
	my ($fn, $de, $colKey, $colValue) = @_;			# $fn : name of file 
	
	open(TOHASH, $fn) or die "Can't open $fn: $!\n";
	
	$_ = <TOHASH>;		# ignore title line
	
	my %retHash = ();
	
	if ($de eq ' ') {		# use special function of "split" by ' '
		while(<TOHASH>) {
			chomp;
			my @d = split ' ';
			my ($k, $v) = @d[$colKey, $colValue];
			
			$retHash{lc $k} = $v;		# lower case
		}
	}else {
		while(<TOHASH>) {
			chomp;
			my @d = split $de;
			my ($k, $v) = @d[$colKey, $colValue];
			
			$retHash{lc $k} = $v;		# lower case
		}
	}
	close TOHASH;

	return \%retHash;
}


########################################
### parse SNP file and create %rLDUp ###
########################################
# return : (refHash, refArray)
#	%{refHash}  = ( rs# of SNP in downstream => rs# of SNP in upstream => undef, )
#	@{refArray} = (	\@(SNP rs num, geneIDs, distance from Genes), ) 
#                 from SNP file
#####
# global variable used : $r2threshold

sub parse_SNPfile_readLD {
	my ($fn) = @_;			# $fn : name of SNP file
	
	my %cLD = ();			# cluster of LD ( rs# => hash of rs# of SNPs in LD )
	             			#		many rs# in a cluster will share same hash of the cluster
	my @snpList = ();		# return variable @{refArray}

	open(SNP, $fn) or die "Can't open $fn: $!";
	$_ = <SNP>;
	
	while(<SNP>) {		# read SNPs in a files based on HapMap project and UCSC files
		chomp;
		my ($thisSNP, undef, $pos, $gidSt, $dGSt, undef, $dEx5St, $dEx3St, @ldData) = split /\t/;
			# $d... means 'distance from'. ex) $dGs : distance from Genes

		#rs#		chr	pos		GeneIDs			from_Genes	from_CDSs	from_exon_5'	from_3'	SNPs_in_LD r2
		#rs3766191	1	1007450	54991			0			686			-780			-389	rs10907177 0.908	rs10907178 0.89
		#rs2272756	1	871896	148398,26155	2072,0		2500,0		,108			,1478	

		my @thisSnpData = ($thisSNP, $gidSt, $dGSt);
		push @snpList, \@thisSnpData;
		
		##
		## - Select SNPs in LD based on R2 value and save the link in hash %cLD
		##
		if(exists $cLD{$thisSNP}) {
			$cLD{$thisSNP}->{$thisSNP} = $pos;
		}
		foreach my $ldDatum (@ldData) {
			chomp $ldDatum;
			my ($snpLDDn, $r2) = split(' ', $ldDatum);
			
			##
			## classification of significant LD by r2 value
			next unless $r2 >= $r2threshold;
			
			## create %cLD
			if(! exists $cLD{$thisSNP} && ! exists $cLD{$snpLDDn}) {
				my %cluster = (
					$thisSNP => $pos,
					$snpLDDn => undef,
				);
				$cLD{$thisSNP} = \%cluster;
				$cLD{$snpLDDn} = \%cluster;		# share hash reference
			}elsif(exists $cLD{$thisSNP} && ! exists $cLD{$snpLDDn}) {
				$cLD{$thisSNP}->{$snpLDDn} = undef;
				$cLD{$snpLDDn} = $cLD{$thisSNP};		# share hash reference
			}elsif(! exists $cLD{$thisSNP} && exists $cLD{$snpLDDn}) {
				$cLD{$snpLDDn}->{$thisSNP} = $pos;
				$cLD{$snpLDDn}->{$snpLDDn} = undef;
				$cLD{$thisSNP} = $cLD{$snpLDDn};		# share hash reference
			}else {
				mergeLdSet (\%cLD, $thisSNP, $snpLDDn);	# share hash reference
			}
		}
	}
	close SNP;
	
	##
	## - Create %rLDUp from %cLD
	##
	my %ldUp = ();

	foreach my $rThisSnpData (@snpList) {
		my ($thisSNP) = @$rThisSnpData;
		
		next unless exists $cLD{$thisSNP};
		
		my $thisPos = $cLD{$thisSNP}->{$thisSNP};
		
		foreach my $ldSnp (keys %{$cLD{$thisSNP}}) {
			# Because there are some missing HapMap SNPs in UCSC data files
			# position of snp in LD couldn't be identified and 
			# the code "defined $cLD{$thisSNP}->{$ldSnp}" has to be in.
			if(defined $cLD{$thisSNP}->{$ldSnp} && $cLD{$thisSNP}->{$ldSnp} < $thisPos) {
				$ldUp{$thisSNP}{$ldSnp} = undef;
			}
		}
	}
	
	return (\%ldUp, \@snpList);

	###
	### merge two exclusive (ex. 'pG') sets of snps in ld
	###
	
	sub mergeLdSet {
		my ($rcld, $snpUp, $snpDn) = @_;
		
		return undef if($rcld->{$snpUp} eq $rcld->{$snpDn});	# identical hash reference
		
		my ($rsrc, $rtrg) = ();
		if( scalar(keys %{$rcld->{$snpUp}}) < scalar(keys %{$rcld->{$snpDn}}) ) {
			($rsrc, $rtrg) = ($rcld->{$snpUp}, $rcld->{$snpDn});
		}else {
			($rsrc, $rtrg) = ($rcld->{$snpDn}, $rcld->{$snpUp});
		}
		
		if(exists $rsrc->{'pG'}) {
			foreach my $exonKey (keys %{delete $rsrc->{'pG'}}) {
				$rtrg->{'pG'}->{$exonKey} = undef;
			}
		}
		
		foreach my $eleSource (keys %$rsrc) {
			$rtrg->{$eleSource} = $rsrc->{$eleSource};		# copy element
			$rcld->{$eleSource} = $rtrg;					# copy hash
		}
		
		%{$rsrc} = ();
		$rsrc = $rtrg;
	}
}

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
