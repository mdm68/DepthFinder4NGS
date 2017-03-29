#!/opt/nasapps/applibs/perl-5.16.2/bin/perl

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use FileHandle;

#Creating the ENSEMBL registry for connecting to ENSEMBL
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
   -host => 'ensembldb.ensembl.org',
   # -host => 'useastdb.ensembl.org',
   -user => 'anonymous'
);


my $filename = $ARGV[0];
open(my $filehandle, '<', $filename) or die "Could not open $filename\n";

my @resultarray;
while(my $line = <$filehandle>){
    chomp $line;
   	@resultarray = split "," , $line;
   	my $fh = FileHandle->new;
    if ($fh->open(">g1.csv")) {
   foreach (@resultarray){
   #file handle 
   		my $gene_adaptor =Bio::EnsEMBL::Registry->get_adaptor( "human", "core","gene" );
		#Inputting the HGNC symbol if interest to obtain the GTE data
		my @genes = @{$gene_adaptor->fetch_all_by_external_name($_)};
		foreach(@genes){
		my $gene=shift(@genes);
    		my $gstring = feature2string($gene);
			#Regex to extract the exact
			foreach($gstring){
				s/\(|\(\(-|\+\d\)//g|s/\(|\(\(-|\+\d//g|s/\)//g|s/\:\s/Chr/g|s/\s//g; # do the replacement
			}	
  				my $transcripts = $gene->get_all_Transcripts();
    					foreach my $transcript (@{$gene->get_all_Transcripts}) {
        					my $tstring = feature2string($transcript);
						foreach($tstring){
						s/\(|\(\(-|\+\d\)//g|s/\(|\(\(-|\+\d//g|s/\)//g|s/\:\s/Chr/g|s/\s//g; # do the replacement
						}
						print $fh "$gstring,$tstring,";
	     				foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
            				my $estring = feature2string($exon);
							foreach($estring){
        						s/\(|\(\(-|\+\d\)//g|s/\(|\(\(-|\+\d//g|s/\)//g|s/\:\s/Chr/g|s/\s//g; # do the replacement
			}
			print $fh "$estring,";
			#written to eliminate the last comma for the last exon, but dint work.
								foreach($estring){
                        			s/\,\\n/\\n/g
                        }
        	}
		print $fh "\n";
          }
    }
print $fh "$_\n\n";
}
}
$fh->close;
}		
exit;
#This subroutine gets and sets the features
    sub feature2string
{
    my $feature = shift;
    my $stable_id  = $feature->stable_id();
    my $seq_region = $feature->slice->seq_region_name();
    my $start      = $feature->start();
    my $end        = $feature->end();
    my $strand     = $feature->strand();

    return sprintf( "%s: %s:%d-%d (%+d)",
        $stable_id, $seq_region, $start, $end, $strand );
}
