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
@genes = @{ $gene_adaptor->fetch_all_by_external_name('ST5') };
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
my $transcripts = $gene->get_all_Transcripts();
while ( my $transcript = shift @{$transcripts} ) {
 

 my @exons = @{$transcript->get_all_Exons()};
 print "\n transcript stable id", $transcript->stable_id, "\n";
 foreach my $exon (@exons) {
  print "\n===================exon start==================================\n";
  #print "\nid: ", $exon->id, ", dbId: ", $exon->dbID, ", display_id: ", $exon->display_id, "\n";
  print "start: ", $exon->start, ", end: ", $exon->end, ", length: ", $exon->length, ",   strand", $exon->strand, "\n";
  print "seq start: ", $exon->seq_region_start, ", seq end: ", $exon->seq_region_end, ", seq length: ", $exon->seq_region_length, ",   seq strand", $exon->seq_region_strand, "\n";
  print "seq: ", $exon->seq->seq, "\n";  
  print "\n====================exon end=================================\n";  
  
  insert_exon($transcript,$host,$exon);
 
 }


}




sub urlencode {
    my $s = shift;
    $s =~ s/ /+/g;
    $s =~ s/([^A-Za-z0-9\+-\/])/sprintf("%%%02X", ord($1))/seg;
    return $s;
}




sub insert_exon{
  $transcript = $_[0];
  $host = $_[1];
  $exon = $_[2];
  $path = $host . "/transcript_exons";
  
  $params = "exon[transcript_stable_id]=" . $transcript->stable_id;
  $params = $params . "&" . exon2string($exon);
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

sub exon2string{

    my $exon = shift;

    my $seq_region_start = $exon->seq_region_start;
    my $seq_region_end = $exon->seq_region_end;
    my $seq_region_strand = $exon->seq_region_strand;
    my $length = $exon->length;
    my $seq = $exon->seq->seq;
    
    my $ret_string = "";
    
    $ret_string = sprintf("exon[seq_region_start]=%s&exon[seq_region_end]=%s&exon[length]=%s&exon[seq_region_strand]=%s&exon[seq]=%s", $seq_region_start, $seq_region_end, $length,$seq_region_strand,$seq);
    return $ret_string;

}