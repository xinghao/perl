use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

#$registry->version_check();

my $host = "http://127.0.0.1:3000";
my $gene_name = 'MEN1';



$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
$va_adaptor = $registry->get_adaptor("human","variation","variationannotation");

@genes = @{ $gene_adaptor->fetch_all_by_external_name($gene_name) };
$gene = @genes[0];




 @vas = @{$va_adaptor->fetch_all_by_associated_gene($gene_name)};
 $count = 0;
 $error_count = 0;
 foreach $va (@vas) {
 	  $count = $count + 1;
    print "Variant Name: ", $va->variation_names(), ", ", $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(),"\n";
    $var = $va->variation();

    print "Class: ", $var->var_class, "\n";
    print "Ambig Code: ", $var->ambig_code, "\n";
    print "Ancestral Code: ", $var->ancestral_allele(), "\n";
    if (associated_gene($var, 'MEN1') > 0) {
      print "Match!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    } else {
      print "Not Match!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    } 
    
    
    
    
    $url = "";
    $url = $url . "variant[variant_type]=dna&";
    $url = $url . "variant[ensembl_id]=". $var->dbID . "&";    
    $url = $url . "variant[variant_class]=". $var->var_class . "&";
    $url = $url . "variant[name]=". $var->name . "&";
    $url = $url . "variant[ambig_code]=". $var->ambig_code . "&";
    $url = $url . "variant[ancestral_allele]=". $var->ancestral_allele . "&";
    $url = $url . "variant[phenotype_name]=". $va->phenotype_name . "&";
    $url = $url . "variant[phenotype_description]=". $va->phenotype_description . "&";
        
    $url = $url . "variant[is_somatic]=". $var->is_somatic . "&";
    $url = $url . "variant[five_prime_flanking_seq]=". $var->five_prime_flanking_seq() . "&";
    $url = $url . "variant[three_prime_flanking_seq]=". $var->three_prime_flanking_seq() . "&";
    
    
    $v_count = 0;
	 foreach $vf (@{$var->get_all_VariationFeatures()}) {
	   $v_count = $v_count + 1;
      print $vf->seq_region_name(), $vf->seq_region_start(), '-',$vf->seq_region_end(), ":", $vf->seq_region_strand, "\n";
      print "region_length: ", $vf->seq_region_length, "\n";
      print "is_somatic: ", $vf->is_somatic, "\n";
      print "display_id: ", $vf->display_id, "\n";
      
  		print "class_SO_term:", $vf->class_SO_term, ", display_consequence: ", $vf->display_consequence, "\n";
  		if($vcount > 1){
  		  print "Duplicated features \n";
  		  $error_count = $error_count + 1;
  		}
  		
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
                     		  		
	  }
	  
	  
    $url =~ s/ /+/g;
    $url =~ s/:/%3a/g;
    $url =~ s/;/+/g;
    
    $path = $host . "/genes/" . $gene->stable_id . "/variants";
    #printf ("URL: " . $path . "/create?". $url);
	  
                 
    insert_variant($url, $path . "?" . $url);  
    

	  
	  last;
    
  }
  
  print "Total: ", $count, ", Error: " , $error_count;

  
sub fetch_variant_by_feature{
  $vf = $_[0];
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
  my @annotations = @{$var->get_all_VariationAnnotations()}
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
  
  
  $url =~ s/ /+/g;
  $url =~ s/:/%3a/g;
  $url =~ s/;/+/g;
  
  $path = $host . "/genes/" . $gene->stable_id . "/variants";
  #printf ("URL: " . $path . "/create?". $url);
  
               
  insert_variant($url, $path . "?" . $url);      
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
  