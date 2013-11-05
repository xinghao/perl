use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $host = "http://127.0.0.1:3000";
$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
$va_adaptor = $registry->get_adaptor("human","variation","variationannotation");
my $gene_name = 'MEN1';
@genes = @{ $gene_adaptor->fetch_all_by_external_name($gene_name) };
$gene = @genes[0];



my $va_adaptor = $registry->get_adaptor('human', 'variation', 'variation'); #get the different adaptors for the different objects needed
my $vf_adaptor = $registry->get_adaptor('human', 'variation', 'variationfeature');



my @rsIds = qw(rs117705251);#COSM23045);#COSM22596);
foreach my $id (@rsIds){
# get Variation object
  my $var = $va_adaptor->fetch_by_name($id); #get the Variation from the database using the name
  get_VariationFeatures($var);
}

sub get_VariationFeatures{
  my $var = shift;
  # get all VariationFeature objects: might be more than 1 !!!
  foreach my $vf (@{$vf_adaptor->fetch_all_by_Variation($var)}){
      print "RSID: ", $vf->variation_name(),"\n"; # print rsID
      print "Allele: ", $vf->allele_string(),"\n"; # print alleles
      print "Consequences: ", join("\n",@{$vf->consequence_type()}),"\n"; # print consequenceType
      #print substr($var->five_prime_flanking_seq,-10) , "[",$vf->allele_string,"]"; #print the allele string
      #print substr($var->three_prime_flanking_seq,0,10), "\n"; # print RefSeq
      print "Sequence name: ", $vf->seq_region_name, ":", $vf->start,"-",$vf->end; # print position in Ref in format Chr:start-end
      #print "gtv return: ", get_TranscriptVariations($vf, $gene); "\n"# get Transcript information
      get_TranscriptVariations($var, $vf, $gene);
  }
}

sub get_hgvs{
  my $tv = shift;
        my $hgvs_transcript = $tv->hgvs_transcript();
        foreach my $key (keys %{$hgvs_transcript}){
          print "hgvs_transcript      ", $key, "->", $hgvs_transcript->{$key}, "\n";
        }

        print "\n\n";
        my $hgvs_protein = $tv->hgvs_protein();
        foreach my $key (keys %{$hgvs_protein}){
          print "hgvs_protein     ", $key, "->", $hgvs_protein->{$key}, "\n";
        }
        
        print "\n\n";
        
        my $hgvs_genomic = $tv->hgvs_genomic();
        foreach my $key (keys %{$hgvs_genomic}){
          print "hgvs_genomic       ", $key, "->", $hgvs_genomic->{$key}, "\n";
        }
  
  
}

sub get_TranscriptVariations1{
  my $vf = shift; 
  
  # get all TranscriptVariation objects: might be more than 1 !!!
  my $transcript_variations = $vf->get_all_TranscriptVariations; #get ALL the effects of the variation in 
                                                                    # different Transcripts
  if (defined $transcript_variations){
    foreach my $tv (@{$transcript_variations}){
      my $gene = $gene_adaptor->fetch_by_transcript_id($tv->transcript->dbID);
      if($gene->stable_id eq "ENSG00000133895") {
      
        print "\n======================================================";
        print "\nPep allele String: ", $tv->pep_allele_string, "\n" if (defined $tv->pep_allele_string);
                                                # the AA change, but only if it is in a coding region
        print "cdna allele string: ", $tv->cdna_allele_string, ", start: ", $tv->cdna_start, ", end: ", $tv->cdna_end, "\n" if (defined $tv->cdna_allele_string);                                              
        print "cds start: ", $tv->cds_start, ", end: ", $tv->cds_end, "\n";
        print "codon_position: ", $tv->codon_position, ", string: ", $tv->codons, "\n";
        get_hgvs($tv);
                                       
        
        print "\n##############################################";
        print "consequence type: ", (join ",", @{$tv->consequence_type}), "\n";
    	print "cdna coords: ", $tv->cdna_start, '-', $tv->cdna_end, "\n";
    	print "cds coords: ", $tv->cds_start, '-', $tv->cds_end, "\n";
    	print "pep coords: ", $tv->translation_start, '-',$tv->translation_end, "\n";
    	print "amino acid change: ", $tv->pep_allele_string, "\n";
    	print "codon change: ", $tv->codons, "\n";
    	print "allele sequences: ", (join ",", map { $_->variation_feature_seq } @{ $tv->get_all_TranscriptVariationAlleles }), "\n";
        print "	distance_to_transcript: ", $tv->distance_to_transcript, "\n";
        print "======================================================\n";
        
        
      }
    }
  }
  print "\n";
}



sub get_TranscriptVariations{
  my $var = $_[0];
  my $vf = $_[1];
  my $gene = $_[2];
  
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
		$url = $host . "/transcript_variants";
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



