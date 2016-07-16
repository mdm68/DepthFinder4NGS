#!/opt/nasapps/applibs/perl-5.16.2/bin/perl

use warnings;
use strict;
use Math::NumberCruncher;
use FileHandle;
my @fields=();

open(DATA1, "<BRAF.csv");

# Open new file to write
#open(DATA2, ">Depths.csv");
 my $fh = FileHandle->new;
 if ($fh->open(">Depths_BRAF.csv")) {

# Copy data from one file to another.
while(my $line = <DATA1>)
{
   #Copying the data from the CSV to new CSV file
   print DATA2 $line;
   @fields = split "," , $line;
   splice @fields,0,2;
   foreach (@fields){
   	         print $fh $_;
		 #Splitting the term at ":" for obtaining the start & the ends eventually
                 my @super_exact_fields=split":",$_;
                 #Regular expression to obtain the chromosome name from the file, it matches the number at the ending of the word
                 ($super_exact_fields[0]) = $super_exact_fields[0]=~/.*?([0-9]+)$/;
                 print "chromosome name check: $super_exact_fields[0]\n" ;	
                 my @final_array=split"-",$super_exact_fields[1];
                 print "ucsc_final3_chr$super_exact_fields[0].wig\n";
                 print "end : $final_array[1]\n";
		#Calling the Subroutine to calculate the depths, by passing the wig file of the interested chromosome, extracted from the CSV file
                 my @array= getDepths("goi_chr$super_exact_fields[0].wig",$final_array[0],$final_array[1]);
   
   				 #calculating stats using API
				my $StdDev = Math::NumberCruncher::StandardDeviation(\@array);
				my ($high,$low) = Math::NumberCruncher::Range(\@array);
				my $mean = Math::NumberCruncher::Mean(\@array);
				my $median = Math::NumberCruncher::Median(\@array);
				my $mode = Math::NumberCruncher::Mode(\@array);

				#printing the stats:
				##print $fh "Std Dev:$StdDev\n", "Range : \n High :$high\n", "Mean:$mean\n","Median:$median\n","Mode:$mode\n";
				print $fh "Min $low Max $high Median $median Mode $mode Mean $mean Stdev $StdDev";
				#print $fh "Std Dev:$StdDev\n", "Range : \n High :$high \n Low: $low\n", "Mean:$mean\n","Median:$median\n","Mode:$mode\n"; 

				#print "Exon Ends,\n\n";                
				print $fh ",";
}
				#print "Transcript Ends,\n\n";                
                                print $fh "\n";

}
#Close the Opened file
close( DATA1 );
}
#close the opened file handle
close( $fh );

#This is the sub-routine to get the depths of the given WIG file of the chromosome of interest
sub getDepths{
my $wig = "/is1/projects/ccrifx/scratch/AllUsers/monica/Project/BAMfromBWA/";
$wig .= shift;
my $start = shift;
my $end= shift;
my @depths=();

	#Open the wig file for the processing to start
	open(WIG, $wig);
	#This loop consists of checking the condition for the start & ends obtained from the WIG File , if they are within the range of the start & end of 
	# the ENSEMBL database
        while(<WIG>) {
	    my $line = $_;
            chomp($line);
	    my @columns = split('\t', $line);
    	    ($columns[0]) = $columns[0]=~/((\d)*)$/;
    	    ($columns[1]) = $columns[1]=~/((\d)*)$/;
            if(($columns[0]>=$start)&&($columns[0]<=$end)){
                push(@depths,$columns[1]);
		print "depth check:".$columns[1]."\n";
        }
        }
        return @depths;       
       }
         close(WIG);
exit;