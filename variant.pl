use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);


$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
$va_adaptor = $registry->get_adaptor("human","variation","variationannotation");
$variation_adaptor = $registry->get_adaptor('human', 'variation', 'variation');
$vfa = $registry->get_adaptor("human","variation","variationfeature");

@genes = @{ $gene_adaptor->fetch_all_by_external_name('MEN1') };
$gene = @genes[0];

  




$var = $variation_adaptor->fetch_by_name('COSM22596');
my $variation_class = $var->var_class();


print "Class: ", $variation_class, "\n";
print "Ambig Code: ", $var->ambig_code, "\n";
print "Ancestral Code: ", $var->ancestral_allele(), "\n";
print "clinical_significance: ", $var->clinical_significance(), "\n";
print "display_consequence: ", $var->display_consequence(), "\n";



@all_syns = @{$var->get_all_synonyms()};
foreach $syn (@all_syns) {
  print: "SYN: ", $syn, "\n";
}

@alleles = @{$var->get_all_Alleles()};
foreach $allele (@alleles) {
  print "allele: " , $allele->allele, ", count: ", $allele->count, ",  subsnp: ", $allele->subsnp, ", subsnp_handle: ", $allele->subsnp_handle, "\n" 
}

@genes = @{$var->get_all_Genes()};
foreach $gene (@genes) {
	print "Gene: ", $gene->external_name, "\n"
}

  foreach $va (@{$va_adaptor->fetch_all_by_Variation($var)}) {
    print $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(),"\n";
  }


 foreach $vf (@{$var->get_all_VariationFeatures()}) {
  #fetch_variant_by_feature($vf, $gene);
  }
 
 
 
 
 
 
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
  
  
  $url =~ s/ /+/g;
  $url =~ s/:/%3a/g;
  $url =~ s/;/+/g;
  
  $path = $host . "/genes/" . $gene->stable_id . "/variants";
  #printf ("URL: " . $path . "/create?". $url);
  
               
#  insert_variant($url, $path . "?" . $url);      
}

