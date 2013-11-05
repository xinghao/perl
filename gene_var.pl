use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

#$registry->version_check();

#my $host = "http://127.0.0.1:3000";
my $host = "http://codondex.com";
my $gene_name = 'MEN1';



$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
$va_adaptor = $registry->get_adaptor("human","variation","variationannotation");

@genes = @{ $gene_adaptor->fetch_all_by_external_name($gene_name) };
$gene = @genes[0];


printf ("seq_region_name: %s, length: %d, start:%d, end:%d\n" , $gene->seq_region_name(), $gene->seq_region_length(),$gene->seq_region_start(),$gene->seq_region_end());


printf ("name: %s, dbId: %s\n, stable_id: %s" , $gene->external_name(), $gene->dbID(), $gene->stable_id);
printf ("sequence name: %s\n", $gene->seqname());
$slice = $gene->slice();
printf ("slice_name: %s\n", $slice->name());


my $vf_adaptor = $registry->get_adaptor('Human', 'variation', 'variationfeature'); #get adaptor to VariationFeature object
my $vfs = $vf_adaptor->fetch_all_by_Slice($slice); #return ALL variations defined in $slice

$icount = 0;
foreach my $vf (@{$vfs}){
  if ($vf->seq_region_start >= $gene->seq_region_start && $vf->seq_region_end <= $gene->seq_region_end && associated_gene($vf->variation,$gene->external_name) > 0) {

    $icount = $icount + 1;
    
    fetch_variant_by_feature($vf, $gene);    
    
    #last;      
      
      
  }
}

printf "Total: ", $icount, "\n";



 sub fetch_variant_by_feature{
  $vf = $_[0];
  $gene = $_[1];
  $var = $vf->variation();
  
  
  print "Class: ", $var->var_class, "\n";
  print "Ambig Code: ", $var->ambig_code, "\n";
  print "Ancestral Code: ", $var->ancestral_allele(), "\n";
  
  
  # for var
  $url = "";
  $url = $url . "variant[variant_type]=dna&";
  $url = $url . "variant[ensembl_id]=". $var->dbID . "&";    
  $url = $url . "variant[variant_class]=". $var->var_class . "&";
  $url = $url . "variant[name]=". $var->name . "&";
  $url = $url . "variant[ambig_code]=". $var->ambig_code . "&";
  $url = $url . "variant[ancestral_allele]=". $var->ancestral_allele . "&";      
  $url = $url . "variant[is_somatic]=". $var->is_somatic . "&";
  $url = $url . "variant[five_prime_flanking_seq]=". $var->five_prime_flanking_seq() . "&";
  $url = $url . "variant[three_prime_flanking_seq]=". $var->three_prime_flanking_seq() . "&";
  
  

  # for va
  my @annotations = @{$var->get_all_VariationAnnotations()};
  foreach $va (@annotations) {
    $url = $url . "variant[phenotype_name]=". $va->phenotype_name . "&";
    $url = $url . "variant[phenotype_description]=". $va->phenotype_description . "&";
    last;    
  }  
  
  # for vfs
  print $vf->seq_region_name(), $vf->seq_region_start(), '-',$vf->seq_region_end(), ":", $vf->seq_region_strand, "\n";
  print "region_length: ", $vf->seq_region_length, "\n";
  print "is_somatic: ", $vf->is_somatic, "\n";
  print "display_id: ", $vf->display_id, "\n";
  print "class_SO_term:", $vf->class_SO_term, ", display_consequence: ", $vf->display_consequence, "\n";
  my $allele_string = $vf->allele_string;
  #$allele_string = ~ s/\// /g;
  print "allele_string: ", $allele_string, "\n";
  
  $url = $url . "variant[allele_string]=" .$allele_string . "&";
  $url = $url . "variant[display_consequence]=" .$vf->display_consequence . "&";
  $url = $url . "variant[class_so_term]=". $vf->class_SO_term . "&";
  $url = $url . "variant[display_id]=". $vf->display_id . "&";
  $url = $url . "variant[length]=". $vf->length . "&";
  $url = $url . "variant[seq_region_end]=". $vf->seq_region_end . "&";
  $url = $url . "variant[seq_region_length]=". $vf->seq_region_length . "&";             
  $url = $url . "variant[seq_region_name]=". $vf->seq_region_name . "&";
  $url = $url . "variant[seq_region_start]=". $vf->seq_region_start . "&";
  $url = $url . "variant[seq_region_strand]=". $vf->seq_region_strand . "&";
  $url = $url . "variant[seqname]=". $vf->seqname . "";
  
  
  my @terms = @{$vf->consequence_type()};
  $i_term = 0;
  foreach $term (@terms) {
    $url = $url . "consequence_types[" . $i_term. "][name]=". $term . "&";
    print "Term: ", $term, "\n";
    $i_term = $i_term + 1;
  }
  
  #add transcript variants
  #$url = $url . get_TranscriptVariations($vf, $gene);
  
  $url =~ s/ /+/g;
  $url =~ s/:/%3a/g;
  $url =~ s/;/+/g;
  
  $path = $host . "/genes/" . $gene->stable_id . "/variants";
  printf ("URL: " . $path . "/create?". $url);
  
               
  insert_variant($url, $path . "?" . $url);
  #add transcript variants
  $path = $host . "/transcript_variants";
  get_TranscriptVariations($var, $vf, $gene, $path);
}


sub insert_variant{
  $params = $_[0];
  $url = $_[1];
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST',$url;
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Error: " . $res->status_line . "\n" }
}



sub associated_gene {
  my $var = $_[0];
  my $gene_name = $_[1];
  
  @genes = @{$var->get_all_Genes()};
  foreach $gene (@genes) {
    printf $gene->external_name , "\n";
    if($gene->external_name eq $gene_name) {
      return 1;
    } 
  } 
  return 0;  
}
  
sub get_TranscriptVariations{
  my $var = $_[0];
  my $vf = $_[1];
  my $gene = $_[2];
  my $url = $_[3];
  
  my $ret_string = "";
  my $i_count = 0;
  
  # get all TranscriptVariation objects: might be more than 1 !!!
  my $transcript_variations = $vf->get_all_TranscriptVariations; #get ALL the effects of the variation in 
                                                                    # different Transcripts
  if (defined $transcript_variations){
    foreach my $tv (@{$transcript_variations}){
      my $gene_t = $gene_adaptor->fetch_by_transcript_id($tv->transcript->dbID);
      print "\ngene-stable: ", $gene->stable_id, "    gene_t->stable_id-id: ", $gene_t->stable_id, "\n";
      if($gene->stable_id eq $gene_t->stable_id) {
      	print "\n=========match============\n";
		my $params = build_transcript_variants_params($var, $tv, $gene);
		
		insert_transcript_variants($params,$url);
      	$icount = $icount + 1;
      	
#        print "\n======================================================";
#        print "\nPep allele String: ", $tv->pep_allele_string, "\n" if (defined $tv->pep_allele_string);
                                                # the AA change, but only if it is in a coding region
#        print "cdna allele string: ", $tv->cdna_allele_string, ", start: ", $tv->cdna_start, ", end: ", $tv->cdna_end, "\n" if (defined $tv->cdna_allele_string);                                              
#        print "cds start: ", $tv->cds_start, ", end: ", $tv->cds_end, "\n";
#        print "codon_position: ", $tv->codon_position, ", string: ", $tv->codons, "\n";
#        get_hgvs($tv);
                                       
        
#        print "\n##############################################";
#        print "consequence type: ", (join ",", @{$tv->consequence_type}), "\n";
#      print "cdna coords: ", $tv->cdna_start, '-', $tv->cdna_end, "\n";
#      print "cds coords: ", $tv->cds_start, '-', $tv->cds_end, "\n";
#      print "pep coords: ", $tv->translation_start, '-',$tv->translation_end, "\n";
#      print "amino acid change: ", $tv->pep_allele_string, "\n";
#      print "codon change: ", $tv->codons, "\n";
#      print "allele sequences: ", (join ",", map { $_->variation_feature_seq } @{ $tv->get_all_TranscriptVariationAlleles }), "\n";
        
#        print "======================================================\n";
        
        
      }
    }
  }
}

sub build_transcript_variants_params {
  my $var = $_[0];
  my $tv = $_[1];
  my $gene = $_[2];

$ret_string = "transcript_variant[variant_stable_id]=" . $var->dbID; 
$ret_string = $ret_string . sprintf("&transcript_variant[stable_id]=%s&transcript_variant[transcript_stable_id]=%s", $tv->dbID, $tv->transcript_stable_id);
$ret_string = $ret_string . sprintf("&transcript_variant[consequence_type]=%s", (join ",", @{$tv->consequence_type}));
$ret_string = $ret_string . sprintf("&transcript_variant[distance_to_transcript]=%s", $tv->distance_to_transcript);
$ret_string = $ret_string . sprintf("&transcript_variant[cdna_start]=%s&transcript_variant[cdna_end]=%s&transcript_variant[cdna_allele_string]=%s", $tv->cdna_start,  $tv->cdna_end, $tv->cdna_allele_string);
$ret_string = $ret_string . sprintf("&transcript_variant[cds_start]=%s&transcript_variant[cds_end]=%s", $tv->cds_start,  $tv->cds_end);
$ret_string = $ret_string . sprintf("&transcript_variant[pep_start]=%s&transcript_variant[pep_end]=%s&transcript_variant[pep_allele_string]=%s", $tv->translation_start,  $tv->translation_end, $tv->pep_allele_string);
$ret_string = $ret_string . sprintf("&transcript_variant[codon_position]=%s&transcript_variant[codons]=%s", $tv->codon_position,  $tv->codons);
$ret_string = $ret_string . sprintf("&transcript_variant[allele_sequences]=%s", (join ",", map { $_->variation_feature_seq } @{ $tv->get_all_TranscriptVariationAlleles }));
return $ret_string;
}


sub insert_transcript_variants {
  $params = $_[0];
  $url = $_[1];
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST',$url;
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Error: " . $res->status_line . "\n" }

}



