use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

#$registry->version_check();





$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
#my $sa = $registry->get_adaptor( 'Human', 'Core', 'Sequence' );
#@genes = @{$gene_adaptor->fetch_all()};
@genes = @{ $gene_adaptor->fetch_all_by_external_name('MEN1') };
$gene = @genes[0];

printf ("seq_region_name: %s, length: %d, start:%d, end:%d\n" , $gene->seq_region_name(), $gene->seq_region_length(),$gene->seq_region_start(),$gene->seq_region_end());

#foreach my $gene (@genes) {
#  printf ("name: %s, dbId: %s\n" , $gene->external_name(), $gene->dbID());
#}

#$gene = $gene_adaptor->fetch_by_dbID(409309);

printf ("name: %s, dbId: %s\n, stable_id: %s" , $gene->external_name(), $gene->dbID(), $gene->stable_id);
printf ("sequence name: %s\n", $gene->seqname());
#printf ("sequence: %s\n", $gene->seq());
$slice = $gene->slice();
printf ("slice_name: %s\n", $slice->name());

#my $dna =
#    ${ $sa->fetch_by_Slice_start_end_strand( $slice, $gene->seq_region_start(), $gene->seq_region_end(), $gene->seq_region_strand ) };
#printf ("DNA: %s\n", $dna);


my $vf_adaptor = $registry->get_adaptor('human', 'variation', 'variationfeature'); #get adaptor to VariationFeature object

# Get a VariationSetAdaptor on the human variation database
#my $vs_adaptor = $registry->get_adaptor('human','variation','variationset');
my $variation_adaptor = $registry->get_adaptor('human', 'variation', 'variation');

# Get all top-level variation sets
#my $top_vss = $vs_adaptor->fetch_all_top_VariationSets();

# Loop over the top level variation sets and recursively print the subsets
foreach my $top_vs (@{$top_vss}) {
  print_set($top_vs);
}




my $slice_iter = $vf_adaptor->fetch_Iterator_by_Slice($slice);

my $variant_count = 0;
LINE: while(my $vf = $slice_iter->next){
            my $variation = $variation_adaptor->fetch_by_name($vf->variation_name);
			@n_genes = @{$vf->get_nearest_Gene()};
			foreach my $n_gene (@n_genes) {
				print "Nearest Gene: ", $n_gene->external_name, " db ID:", $n_gene->stable_id, "\n";
				if ($n_gene->stable_id ne "ENSG00000133895") {
					print "skip\n";
				
				} else {
		              print "Class: ", $vf->var_class, " Variation: ", $vf->variation_name, " with alleles ", $vf->allele_string, 
		                " in chromosome ", $slice->seq_region_name, " and position ", $vf->start,"-",$vf->end,"\n";
		              
					
				    @synonyms = @{$vf->get_overlapping_Genes()};
				    print "Genes: @synonyms\n";
		    
		              # get Variation object
		              
		              get_VariationFeatures($variation);
		              
		                $variant_count = $variant_count + 1;
		                if ($variant_count > 1) {
		                  last;
		                }
				
				
				}
			}

              
}
               
#my $vfs = $vf_adaptor->fetch_all_by_Slice($slice); #return ALL variations defined in $slice

#foreach my $vf (@{$vfs}){
#                $variant_count = $variant_count + 1;
#                if ($variant_count > 10) {
#                  break;
#                }
#  print "Variation: ", $vf->variation_name, " with alleles ", $vf->allele_string, 
#        " in chromosome ", $slice->seq_region_name, " and position ", $vf->start,"-",$vf->end,"\n";
#}

#@dbentries = @{ $gene->get_all_DBEntries() };
#foreach my $dbentry (@dbentries) {
#  printf("description: %s\n", $dbentry->description())
#}


#printf ("%s\n", $gene->seq());
#printf "hello" . @gene->description("fFF")





# We define a function that will help us recurse over the set hierarchy and print the data   
sub print_set {
  my $set = shift;
  my $indent = shift || "";
  
  # Print the set attributes
  printf("\%s\%s (\%s)\n",$indent,$set->name(),$set->short_name());
  
  # Get the subsets that have the current set as immediate parent
  my $subsets = $set->get_all_sub_VariationSets(1);
  
  # Call the print subroutine for each of the subsets with an increased indentation
  foreach my $subset (@{$subsets}) {
    print_set($subset,$indent . "  ");
  }
}


sub get_VariationFeatures{
  my $var = shift;
  # get all VariationFeature objects: might be more than 1 !!!
  foreach my $vf (@{$vf_adaptor->fetch_all_by_Variation($var)}){
      print $vf->variation_name(),","; # print rsID
      print $vf->allele_string(),","; # print alleles
      print join(",",@{$vf->consequence_type()}),","; # print consequenceType
      print substr($var->five_prime_flanking_seq,-10) , "[",$vf->allele_string,"]"; #print the allele string
      print substr($var->three_prime_flanking_seq,0,10), ","; # print RefSeq
      print $vf->seq_region_name, ":", $vf->start,"-",$vf->end; # print position in Ref in format Chr:start-end
      get_TranscriptVariations($vf); # get Transcript information
  }
}

sub get_TranscriptVariations{
  my $vf = shift; 
  
  # get all TranscriptVariation objects: might be more than 1 !!!
  my $transcript_variations = $vf->get_all_TranscriptVariations; #get ALL the effects of the variation in 
                                                                    # different Transcripts
  if (defined $transcript_variations){
    foreach my $tv (@{$transcript_variations}){
      print ",", $tv->pep_allele_string if (defined $tv->pep_allele_string);
                                              # the AA change, but only if it is in a coding region
      my $gene = $gene_adaptor->fetch_by_transcript_id($tv->transcript->dbID);
      print ",",$gene->stable_id if (defined $gene->external_name); # and the external gene name
    }
  }
  print "\n";
}

