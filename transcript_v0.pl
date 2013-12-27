use Bio::EnsEMBL::Registry;
use LWP::UserAgent;


my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);




#my $host = "http://127.0.0.1:3000";
my $host = "http://www.codondex.com";


$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
@genes = @{ $gene_adaptor->fetch_all_by_external_name('ZNF692') };
$gene = @genes[0];

printf ("name: %s, dbId: %s\n, stable_id: %s" , $gene->external_name(), $gene->dbID(), $gene->stable_id);

printf ("Biutype: %s\n", $gene->biotype());
printf ("dbId: %s\n", $gene->dbID());
printf ("Description: %s\n", $gene->description());
printf ("desplay ID: %s\n", $gene->display_id());
printf ("External Name: %s\n", $gene->external_name());
printf ("Status: %s\n", $gene->status());
printf ("Source: %s\n", $gene->source());
printf ("Species: %s\n", $gene->species());
 printf("Stable_id: %s\n", $gene->stable_id());
printf("Start: %s, End: %s, Length: %s, Strand: %s", $gene->start, $gene->end, $gene->length, $gene->strand);
printf ("seq_name: %s \n" , $gene->seqname());
printf ("seq_region_name: %s, length: %d, start:%d, end:%d\n" , $gene->seq_region_name(), $gene->seq_region_length(),$gene->seq_region_start(),$gene->seq_region_end());




# transcript
my $icount = 0;
@ts = @{$gene->get_all_Transcripts()};
my $t_size = @ts;
printf("\n Size: %s \n", $t_size);
  
my $transcripts = $gene->get_all_Transcripts();

while ( my $transcript = shift @{$transcripts} ) {
$icount = $icount + 1;
 insert_transcript($transcript,$host,$gene); #transcript2string($transcript);


}

printf("\n Size: %s \n", $icount);



sub insert_gene{
  $params = $_[0];
  $path = $_[1];
  
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST',$path;
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Error: " . $res->status_line . "\n" }
}

sub urlencode {
    my $s = shift;
    $s =~ s/ /+/g;
    $s =~ s/([^A-Za-z0-9\+-\/])/sprintf("%%%02X", ord($1))/seg;
    return $s;
}

sub insert_transcript{
  $transcript = $_[0];
  $host = $_[1];
  $gene = $_[2];
  
  $path = $host . "/transcripts";
  
  $params = "transcript[gene_stable_id]=" . $gene->stable_id;
  $params = $params . "&" . transcript2string($transcript);
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST',$path;
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Transcript Error: " . $res->status_line . "\n" }


}



sub transcript2string
{
    my $transcript = shift;

    my $stable_id  = $transcript->stable_id();
    my $seq_region = $transcript->slice->seq_region_name();
    my $start      = $transcript->start();
    my $end        = $transcript->end();
    my $strand     = $transcript->strand();
    
    my $ret_string = "";
    
    printf( "\n%s: %s:%d-%d (%+d)",
        $stable_id, $seq_region, $start, $end, $strand );
    $ret_string = sprintf("transcript[stable_id]=%s&transcript[seq_region_name]=%s&transcript[seq_region_start]=%s&transcript[seq_region_end]=%s&transcript[seq_region_length]=%s&transcript[seq_region_strand]=%s&transcript[seq]=%s",$stable_id, $seq_region, $start, $end, $transcript->seq_region_length,$strand,$transcript->seq->seq);
    #printf("\nseq_region_start: %s, end: %s, strand: %s", $transcript->seq_region_start, $transcript->seq_region_end,$transcript->seq_region_strand);
    
    printf("\ncDna_start: %s, end: %s, seq: %s", $transcript->cdna_coding_start, $transcript->cdna_coding_end,$transcript->spliced_seq);
    $ret_string = $ret_string . "&" . sprintf("transcript[cdna_start]=%s&transcript[cdna_end]=%s&transcript[cdna_seq]=%s",$transcript->cdna_coding_start, $transcript->cdna_coding_end,$transcript->spliced_seq);
    printf("\ncCds_start: %s, end: %s, seq: %s", $transcript->coding_region_start, $transcript->coding_region_end,$transcript->translateable_seq);
    $ret_string = $ret_string . "&" . sprintf("transcript[cds_start]=%s&transcript[cds_end]=%s&transcript[cds_seq]=%s",$transcript->coding_region_start,$transcript->coding_region_end,$transcript->translateable_seq);
    
    my $protein = $transcript->translate();
    if (defined($protein)){
      printf("\ncProtein seq: %s", $protein->seq);
      $ret_string = $ret_string . "&" . sprintf("transcript[protein_seq]=%s",$protein->seq);
    }else{
      printf("\n No Protein Exists")
    }    
    
    
    printf("\n");
    return $ret_string;
      
}




sub insert_intron{
  $transcript = $_[0];
  $host = $_[1];
  $intron = $_[2];
  $path = $host . "/introns";
  
  $params = "intron[transcript_stable_id]=" . $transcript->stable_id;
  $params = $params . "&" . intron2string($intron);
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST',$path;
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Transcript Error: " . $res->status_line . "\n" }

  
}

sub intron2string{

    my $intron = shift;

    my $seq_region_start = $intron->seq_region_start;
    my $seq_region_end = $intron->seq_region_end;
    my $seq_region_strand = $intron->seq_region_strand;
    my $length = $intron->length;
    my $seq = $intron->seq;
    
    my $ret_string = "";
    
    $ret_string = sprintf("intron[seq_region_start]=%s&intron[seq_region_end]=%s&intron[length]=%s&intron[seq_region_strand]=%s&intron[seq]=%s", $seq_region_start, $seq_region_end, $length,$seq_region_strand,$seq);
    return $ret_string;

}