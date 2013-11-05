use URI::Encode;


my $encoder = URI::Encode->new({double_encode => 1});
$url = "http://perl.com/foo bar:ff";
printf urlencode($url), "\n"; # prints http://perl.com/foo%20bar



sub urlencode {
    my $s = shift;
    $s =~ s/ /+/g;
    $s =~ s/([^A-Za-z0-9\+-\/])/sprintf("%%%02X", ord($1))/seg;
    return $s;
}

