use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

#$registry->version_check();





$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
$va_adaptor = $registry->get_adaptor("human","variation","variationannotation");
  
@genes = @{ $gene_adaptor->fetch_all_by_external_name('MEN1') };
$gene = @genes[0];

printf ("seq_region_name: %s, length: %d, start:%d, end:%d\n" , $gene->seq_region_name(), $gene->seq_region_length(),$gene->seq_region_start(),$gene->seq_region_end());


printf ("name: %s, dbId: %s\n, stable_id: %s" , $gene->external_name(), $gene->dbID(), $gene->stable_id);
printf ("sequence name: %s\n", $gene->seqname());
#printf ("sequence: %s\n", $gene->seq());
$slice = $gene->slice();
printf ("slice_name: %s\n", $slice->name());



 @vas = @{$va_adaptor->fetch_all_by_associated_gene('MEN1')};
 $count = 0;
 foreach $va (@vas) {
 	$count = $count + 1;
    print "Variant Name: ", $va->variation_names(), ", ", $va->phenotype_name(), $va->phenotype_description(), $va->source_name(), $va->study_type(),"\n";
    $var = $va->variation();
    
	 foreach $vf (@{$var->get_all_VariationFeatures()}) {
#	    print $vf->seq_region_name(), $vf->seq_region_start(), '-',$vf->seq_region_end(),"\n";
		print "class_SO_term:", $vf->class_SO_term, ", display_consequence: ", $vf->display_consequence, "\n";
#		my @terms = $vf->consequence_type();
#		foreach $term (@terms) {
#			print "Term: ", $term, "\n";
#		} 	
	  }
    
  }
  
  print "Total: ", $count
