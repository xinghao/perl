
use LWP::UserAgent;

insert_gene("match=www&errors=0");

sub insert_gene{
  $params = $_[0];
  
  $ua = LWP::UserAgent->new;
  
  my $req = new HTTP::Request 'POST','http://127.0.0.1:3000/genes/insert';
  $req->content_type('application/x-www-form-urlencoded');
  $req->content($params);

  my $res = $ua->request($req);
        
  
  # check the outcome
  print $res->decoded_content;
  if ($res->is_success) { 
  #print $res->decoded_content; 
  }
  else { print "Error: " . $res->status_line . "\n" }
}