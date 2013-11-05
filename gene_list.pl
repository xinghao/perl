use Bio::EnsEMBL::Registry;
use LWP::UserAgent;
#use URI::Encode;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
);

#$registry->version_check();


#my $host = "http://127.0.0.1:3000";
my $host = "http://codondex.com";


$gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );

@genes = @{$gene_adaptor->fetch_all()};

my $gene_name = "";
my $count = 0;
 foreach $gene (@genes) {
  @ts = @{$gene->get_all_Transcripts()};
  my $t_size = @ts;
  if($t_size > $count) {
    $gene_name = $gene->external_name();
    $count = $t_size;
  }
  #printf ("name: %s, dbId: %s, stable_id: %s, %s\n" , $gene->external_name(), $gene->dbID(), $gene->stable_id, $t_size);
  printf ("%s, %s, %s, %s, %s\n" , $gene->external_name(), $gene->dbID(), $gene->stable_id, $gene->length, $t_size);
  
 }