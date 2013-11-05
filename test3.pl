sub associated_gene {
  my $var = $_[0];
  my $gene_name = $_[1];

  if ($var eq $gene_name) {
    return 1;
  }  
  return 0;  
}
  
  
  
if (associated_gene("abc", "cd1") > 0) {
  printf "match!!!!!!\n";
  }else{
  printf "not match!!!!!!\n";
  } 