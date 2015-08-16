open(fl,"data");
while(<fl>){
	chomp($_);
	@fds = split(/ /, $_);
	if($fds[1]>=0 && scalar(@fds)>3){
		$graph{$fds[0]}{$fds[1]}++;
	}
	
}
close(fl);
open(fl,"author_doc");
while(<fl>){
	chomp($_);
	@fds = split(/\t/, $_);
	if(scalar(@fds)>1){
		for($i=2;$i<scalar(@fds);$i++){
			push(@{$author_map{$fds[0]}}, $fds[$i]);
		}
	}	
}
close(fl);

open(fl,">doc_cited_authors");
foreach $doc (keys %graph){
	foreach $c (keys %{$graph{$doc}}){
		foreach $a (@{$author_map{$c}}){
			push (@{$cited_map{$doc}}, $a) unless grep /$a/, @{$cited_map{$doc}};
		}
	}
}
print fl "1\n";
foreach $doc (keys %cited_map){
	print fl "$doc\t";
	print fl scalar(@{$cited_map{$doc}});
	foreach $a (@{$cited_map{$doc}}){
		print fl "\t$a";
	}
	print fl "\n";
}


