use Lingua::Stem;
$stemmer = Lingua::Stem->new(-locale => 'EN-UK');
$stemmer->stem_caching({ -level => 2 });
my @stopwords=qw(a b c d e f g h i j k l m n o p q r s t u v w x y z about above according across actually adj after afterwards again against all almost alone along already also although always among amongst an and another any anyhow anyone anything anywhere are aren't around as at be became because become becomes becoming been before beforehand begin beginning behind being below beside besides between beyond billion both but by can can't cannot caption co company corp corporation could couldn't did didn't do does doesn't don't down during each eg eight eighty either else elsewhere end ending enough etc even ever every everyone everything everywhere except few fifty first five for former formerly forty found four from further had has hasn't have haven't he he'd he'll he's hence her here here's hereafter hereby herein hereupon hers herself him himself his how however hundred i i'd i'll i'm i've ie if in inc indeed instead into is isn't it it's its itself last later latter latterly least less let let's like likely ltd made make makes many maybe me meantime meanwhile might million miss more moreover most mostly mr mrs much must my myself namely neither never nevertheless next nine ninety no nobody none nonetheless noone nor not nothing now nowhere of off often on once one one's only onto or other others otherwise our ours ourselves out over overall own per perhaps rather recent recently same seem seemed seeming seems seven seventy several she she'd she'll she's should shouldn't since six sixty so some somehow someone something sometime sometimes somewhere still stop such taking ten than that that'll that's that've the their them themselves then thence there there'd there'll there're there's there've thereafter thereby therefore therein thereupon these they they'd they'll they're they've thirty this those though thousand three through throughout thru thus to together too toward towards trillion twenty two under unless unlike unlikely until up upon us used using very via ve was wasn't we we'd we'll we're we've well were weren't what what'll what's what've whatever when whence whenever where where's whereafter whereas whereby wherein whereupon wherever whether which while whither who who'd who'll who's whoever whole whom whomever whose why will with within without won't would wouldn't yeah yes yet you you'd you'll you're you've your yours yourself yourselves 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 http www et al introduction abstract experiments experiment figure reference references acknowledgement acknowledgements conclusion conclusions future work remark thank thanks grant grants use used present presented paper research given result based does finally following good using zhang zha); 

my %stop = map { lc $_ => 1 } @stopwords;

sub findwords {
  my $string = shift;
  my (@ok, %seen);
  while ($string =~ /((\w|')+)/g) {
    push @ok, $1 unless $stop{lc $1};# or $seen{lc $1}++;
  }
  return @ok;
}
if($#ARGV<2){print "USAGE: perl <file> <file format (field number where the text start content_file> <Out file>";}

open(fl,@ARGV[0]);
@lines=<fl>;
close(fl);
open(fl,">@ARGV[2]");
foreach(@lines){
  chomp($_);
  @flds=split(/:/,$_);
  $line="";
  for($i=(@ARGV[1]-1);$i<scalar(@flds);$i++){
    $line=join(" ",$line,$flds[$i]);
  }
  $line =~ s/[^A-Za-z\.:,]/ /g;
  $line =~ s/\s+/ /g;
  $line = join(" ", findwords(lc $line));
  #@wrds=split(/ /,$line);
  #$stemmed_words = $stemmer->stem(findwords(lc $line));
  for($i=0;$i<(@ARGV[1]-1);$i++){
    print fl "$flds[$i]:";
  }
  #foreach(@$stemmed_words){
  print fl "$line\n";
  #}
}
close(fl);

