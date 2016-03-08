use strict;
use warnings;

# ------------------------------------------------------------------------------
# Function: Add GeneID and symbol to refSeq file
# ------------------------------------------------------------------------------
#
# Usage   : >perl preprocess1_GeneIDnSym.pl <Arguments>
#
# INPUT files
#	1: refGene_UCSC_hg18.txt 
#                          http://genome.ucsc.edu/cgi-bin/hgTables?command=start
#                          Table Browser; 
#                             track: RefSeq Genes -> table: refGene
#
#	2: gene2refseq  <= 2008-09-17 ftp://ftp.ncbi.nih.gov/gene/DATA/
#
#	3: gene_info	<= 2008-09-17 ftp://ftp.ncbi.nih.gov/gene/DATA/
#
# OUTPUT file
#	: .\refGeneIDnSym_UCSC_hg18.txt
#
#	0     	1   	2    	3     	4      	5    	6       7    	 	8      	9        	10
#	GeneID	Symbol	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	id	name2	cdsStartStat	cdsEndStat	exonFrames
#	79501	OR4F5	NM_004	chr1	+		58953	59871	58953		59871	1	 		58953,		59871,		0	OR4F5	cmpl	cmpl	0,
#
# ------------------------------------------------------------------------------

my $startTime = time;

sub showUsage {
	print <<END;
Usage: >perl preprocess1_GeneIDnSym.pl (refG) (g2sq) (gIfo) [outF=ref...]

     refG - refGene file name from UCSC table browser
     g2sq - 'gene2refseq' file from ftp://ftp.ncbi.nih.gov/gene/DATA/
     gIfo - 'gene_info' file from ftp://ftp.ncbi.nih.gov/gene/DATA/
     outF - output file name (default = \'.\\refGeneIDnSym_UCSC_hg18.txt\')
END
	exit;
}
showUsage unless ( scalar(@ARGV) >= 3 );

print "Program is running.\n";

## file names
my $inRefFN = shift(@ARGV);
my $inGeneIDFN = shift(@ARGV);
my $inGeneSymbolFN = shift(@ARGV);

my $outFN  = shift(@ARGV) || ".\\refGeneIDnSym_UCSC_hg18.txt";

my %id2sym = ();
{
	open(GNSYM, $inGeneSymbolFN) or die "Can't open $inGeneSymbolFN: $!\n";
	<GNSYM>;
	while(<GNSYM>) {
		chomp;
		my ($tax_id, $geneID, $symbol, undef) = split (/\t/, $_, 4);
		next if ($tax_id != 9606);
		
		#Format: tax_id GeneID Symbol LocusTag Synonyms dbXrefs chromosome map_location description type_of_gene Symbol_from_nomenclature_authority Full_name_from_nomenclature_authority Nomenclature_status Other_designations Modification_date (tab is used as a separator, pound sign - start of a comment)
		#7	5692769	NEWENTRY	-	-	-	-	-	Record to support submission of GeneRIFs for a gene not in Entrez Gene (Azorhizobium caulinodans Dreyfus et al. 1988; Azotirhizobium caulinodans.  Use when strain, subtype, isolate, etc. is unspecified, or when different from all specified ones in Gene.).	other	-	-	-	-	20071023
		
		$id2sym{$geneID} = $symbol;
	}
	close GNSYM;
}

my %rna2id = ();
{
	open(GNID, $inGeneIDFN) or die "Can't open $inGeneIDFN: $!\n";
	<GNID>;
	while(<GNID>) {
		chomp;
		my ($tax_id, $geneID, undef, $rnaNt, undef) = split (/\t/, $_, 5);
		next if ($tax_id != 9606 or $rnaNt eq '-');
		
		#Format: tax_id GeneID status RNA_nucleotide_accession.version RNA_nucleotide_gi protein_accession.version protein_gi genomic_nucleotide_accession.version genomic_nucleotide_gi start_position_on_the_genomic_accession end_position_on_the_genomic_accession orientation assembly (tab is used as a separator, pound sign - start of a comment)
		#9	1246500	PROVISIONAL	-	-	NP_047184.1	10954455	NC_001911.1	10954454	348	1190	-	-
		
		my ($rnaAcc, $rnaVer) = split (/\./, $rnaNt);
		
		$rna2id{$rnaAcc} = $geneID;
	}
	close GNID;
}

{
	open(OUT , ">$outFN") or die "Can't open $outFN: $!\n";
	open(REF, $inRefFN) or die "Can't open $inRefFN: $!\n";
	{
		$_ = <REF>;
		chomp;
		my (undef, @d) = split /\t/;
		
		print OUT join("\t", "GeneID", "Symbol", @d) . "\n";		# append column titles at the beginning
	}
	
	while(<REF>) {
		chomp;
		my (undef, $rna, @d) = split /\t/;
		
		if (exists($rna2id{$rna})) {
			my $id = $rna2id{$rna};
			if (exists($id2sym{$id})) {
				print OUT join("\t", $id, $id2sym{$id}, $rna, @d) . "\n";
			}else {
				print $id . " - no symbol\n";
				print OUT join("\t", $id, '', $rna, @d) . "\n";
			}
		}else {
			print $rna . " - no GeneID\n";
		}
	}
	close REF;
	close OUT;
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