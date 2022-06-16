#!/usr/bin/perl

#process gz seq files

use File::Temp ();
use Archive::Extract ();


@PE1barcodes=qw(
ATCGTACG
);
@PE2barcodes=qw(
CATAGCCG
);
$numPE1barcodes=@PE1barcodes;
$numPE2barcodes=@PE2barcodes;

for($i=0;$i<$numPE1barcodes;$i++) {
	for($j=0;$j<$numPE2barcodes;$j++) {
		$filename="1"."$PE1barcodes[$i]"."2"."$PE2barcodes[$j]";
		#do not use 1st base in barcode because many of them are N, start with 1 instead of 0
		$PE1_7bp=substr($PE1barcodes[$i], 1, 7);
		$PE2_7bp=substr($PE2barcodes[$j], 1, 7);
		$barcodes="$PE1_7bp"."$PE2_7bp";
		$barcodes2filename{$barcodes}=$filename;
	}
}


@files2read = qw(
001
);


$filenametemplate="../../ACH2-fixed-3LTR_S33_L001_";



for my $subnum ( @files2read ) {

	my $R1="$filenametemplate"."R1_"."$subnum".".fastq.gz";
	my $R2="$filenametemplate"."R2_"."$subnum".".fastq.gz";
	
	my $read1 = File::Temp->new();
	my $read2 = File::Temp->new();

  	my $ae1 = Archive::Extract->new( archive => $R1 );
  	print "$R1\n";
  	my $ok1 = $ae1->extract( to => $read1->filename );
  	die $ae1->error unless $ok1;
  	
  	my $ae2 = Archive::Extract->new( archive => $R2 );
  	print "$R2\n";
  	my $ok2 = $ae2->extract( to => $read2->filename );
  	die $ae2->error unless $ok2;

  # do what you want to the uncompressed content in the
  # file named $tmp->filename
  #...
  # at the end of the scope, $tmp vanishes 
  
  
	open INread1, "<$read1" or die;
	open INread2, "<$read2" or die;
	while ($lineread1=<INread1>) {
		chomp($lineread1);
		if($lineread1=~/(^@.*)( 1.*)/) {
		    $read1_fastq_id=$lineread1;
			$totalreads++;
			$read1_id=$1;
			
			$read1_id =~s/\@//;
			$read1_seq=<INread1>;
			chomp($read1_seq);
			<INread1>;  #the + sign in between seq and quality
			$read1_qual=<INread1>;
			chomp($read1_qual);
			
			$lineread2=<INread2>;
			chomp($lineread2);
			$read2_fastq_id = $lineread2;
			($read2_id,$left)=split(/\ /, $lineread2);
			$read2_i7=substr($left,7,7);
			$read2_id =~s/\@//;
			$read2_id_constant = $read2_id;
			$read2_seq=<INread2>;
			chomp($read2_seq);
			<INread2>;
			$read2_qual=<INread2>;
			chomp($read2_qual);
			
			#split read2, 6bp_barcode_linkernestedPrimer+1bpA_junction
			#all linker add a T at the end after ligation
			$check = 0;
			if($read2_seq=~/(..........)(TAGTGC.*TTAGAGGACT)(.*)/) {  #all linker add a T at the end after ligation
				$dogtag = $1;
				$read2_id = $read2_id . "&" . $dogtag;
				$read2_seq_trimmed=$3;
				$read2_qual_start = length($1) + length($2);
				$read2_qual_trimmed = substr($read2_qual, $read2_qual_start);
				$check = 1;
			}
			else{
				$read2_seq_trimmed=$read2_seq;
				$read2_qual_trimmed=$read2_qual;
			}

	
			if($read1_id ne $read2_id_constant){
				print "read1_id does not match read2_id:\n$read1_id\n$read2_id\n";
			}
			
			
			$read1_nestedPrimer=substr($read1_seq,14,23);
			$read1_barcode=substr($read1_seq,7,7);
			$barcodes="$read1_barcode"."$read2_i7";
			$filename=$barcodes2filename{$barcodes};

			$LTR_junction=$read1_seq;
			$read1_qual_trimmed = $read1_qual;
			#if lenti-3'LTR, split like this  nested primer: CCCTTTTAGTCAGTGTGGAAAATC
			#Junction: 3LTRnestedPrimer_7bpLTRend_junction: CCCTTTTAGTCAGTGTGGAAAATC_TCTAGCA_GTAGTAGTT
			if ($read1_nestedPrimer eq "CCTTTTAGTCAGTGTGGAAAATC") { 
				($preLTRend, $postLTRend) = split (/TCTAGCA/, $LTR_junction);
				$LTR_junction= "$preLTRend"."TCTAGCA\t"."$postLTRend";
				$read1_qual_start = length($preLTRend) + 7;
				$read1_qual_trimmed = substr($read1_qual,$read1_qual_start);
			}

			#print format:  id	read1or2	LTR-TCTCTAGCA	insert	qualification
            if ($check) {
			push (@$filename, "$read1_id\t1\t$LTR_junction\t$read1_qual_trimmed\t$read1_fastq_id\t\t\t\t\t\t\t\t$read2_id\t2\t$read2_seq\t$read2_seq_trimmed\t$read2_qual_trimmed\t$read2_fastq_id\t\t\t\t\t\t\t\n");
			}


		}
		
	}
	
	close INread1;
	close INread2;  
}

print "Totalreads:\t$totalreads\n";
while(($barcodes, $filename)=each(%barcodes2filename)) {
	$counts=@$filename;
	print "$filename\:\t$counts\n";
	$file2print="$filename".".txt";
	open OUT, ">$file2print" or die;
	foreach $line (@$filename) {
		print OUT "$line";
	}
	close OUT;
}
	
