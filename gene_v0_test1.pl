use Bio::EnsEMBL::Registry;
use LWP::UserAgent;
#use URI::Encode;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

#$registry->version_check();


my $host = "http://127.0.0.1:3000";
#my $host = "http://www.codondex.com";


$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
#my $sa = $registry->get_adaptor( 'Human', 'Core', 'Sequence' );
#@genes = @{$gene_adaptor->fetch_all()};
@genes = @{ $gene_adaptor->fetch_all_by_external_name('MEN1') };
$gene = @genes[0];

genesize($gene, $host);
printf("\n==================================\n");



sub genesize{
  my $gene =  $_[0];
  my $host =  $_[1];
  printf("\n Gene: %s, size: %s", $gene->external_name(), $gene->length());
  
  # transcript
  my $cdna_max = 0;
  my $cdna_avg = 0;
  my $cds_max = 0;
  my $cds_avg = 0;  
  my $protein_max = 0;
  my $protein_avg = 0;
  my $intron_max = 0;
  my $intron_avg = 0;
  my $exon_max = 0;
  my $exon_avg = 0;          
  my $icount = 0;
  
  
  my $transcripts = $gene->get_all_Transcripts();
  while ( my $transcript = shift @{$transcripts} ) {
   $icount = $icount + 1;
   my $cdna_l = length($transcript->spliced_seq); 
   if ( $cdna_l > $cdna_max ) {$cdna_max = $cdna_l;}
   $cdna_avg = $cdna_avg + $cdna_l;
   my $cds_l = length($transcript->translateable_seq); 
   if ( $cds_l > $cds_max ) {$cds_max = $cds_l;}
   $cds_avg = $cds_avg + $cds_l;     
    my $protein = $transcript->translate();
    if (defined($protein)){
     my $protein_l = length($protein->seq);
     if ( $protein_l > $protein_max ) {$protein_max = $protein_l;}
     $protein_avg = $protein_avg + $protein_l;      
    }    
    
    my @introns = @{$transcript->get_all_Introns()};
    $intron_l = 0;
    foreach my $intron (@introns) {
      $intron_l = $intron_l + $intron->length;            
    }
    if ( $intron_l > $intron_max ) {$intron_max = $intron_l;}
    $intron_avg = $intron_avg + $intron_l; 
    
    my @exons = @{$transcript->get_all_Exons()};
    $exon_l = 0;
    foreach my $exon (@exons) {
      $exon_l = $exon_l + $exon->length;
    }
    if ( $exon_l > $exon_max ) {$exon_max = $exon_l;}
    $exon_avg = $exon_avg + $exon_l; 
    
    
        
  }  
  $cdna_avg = $cdna_avg * 1.0 / $icount;
  $cds_avg = $cds_avg * 1.0 / $icount;
  $protein_avg = $protein_avg * 1.0 / $icount;
  $intron_avg = $intron_avg * 1.0 / $icount;
  $exon_avg = $exon_avg * 1.0 / $icount;
  
  printf("\ncDna[Max: %s, Avg: %s]", $cdna_max, $cdna_avg);
  printf("\nCDS[Max: %s, Avg: %s]", $cds_max, $cds_avg);
  printf("\nProtein[Max: %s, Avg: %s]", $protein_max, $protein_avg);
  printf("\nIntron[Max: %s, Avg: %s]", $intron_max, $intron_avg);
  printf("\nExon[Max: %s, Avg: %s]", $exon_max, $exon_avg);
  
  $path = $host . "/gene_sizes";
  
  
  $params = "gene_size[gene_id]=" . $gene->stable_id;
  $params = $params . "&" . "gene_size[gene_name]=" . $gene->external_name();
  $params = $params . "&" . "gene_size[gene_length]=" . $gene->length();
  $params = $params . "&" . "gene_size[transcripts_amount]=" . $icount;  
  $params = $params . "&" . "gene_size[cdna_max]=" . $cdna_max . "&" . "gene_size[cdna_avg]=" . $cdna_avg;
  $params = $params . "&" . "gene_size[cds_max]=" . $cds_max . "&" . "gene_size[cds_avg]=" . $cds_avg;
  $params = $params . "&" . "gene_size[protein_max]=" . $protein_max . "&" . "gene_size[protein_avg]=" . $protein_avg;
  $params = $params . "&" . "gene_size[intron_max]=" . $intron_max . "&" . "gene_size[intron_avg]=" . $intron_avg;
  $params = $params . "&" . "gene_size[exon_max]=" . $exon_max . "&" . "gene_size[exon_avg]=" . $exon_avg;
  
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST',$path;
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Gene size Error: " . $res->status_line . "\n" }
  
}


