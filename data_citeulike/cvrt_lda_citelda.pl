if($#ARGV<1){print "USAGE:perl <citation_file> <content_file>";}
open(fl,@ARGV[0]);
@citations=<fl>;
close(fl);
open(fl_doc,">docmap");
foreach $citation(@citations){
	chomp($citation);
	@flds=split(/:/,$citation);
	if(($flds[0]=~m/[0-9]/)&&($flds[1]=~m/[0-9]/)){
		$doc{$flds[0]}++;
		$doc{$flds[1]}++;
        	$cited_docs{$flds[0]}++;
	}
}
### Produce document mapping 
$i=0;
%wrd=();
open(fl,@ARGV[1]);
@content=<fl>;
close(fl);
foreach(@content){
        chomp($_);
        @flds=split(/:/,$_);
	if(($flds[0]=~m/[0-9]/)){
		$doc{$flds[0]}++;
        	$line="";
        	for($i=1;$i<scalar(@flds);$i++){
          		$line=join(" ",$line,$flds[$i]);
        	}
        	@wrds=split(/ /,$line);
        	foreach(@wrds){
                	$wrd{$_}++;
        	}
	}
}

$i=0;
%docmap=map{$i++ => $_} (keys %doc);
foreach(keys %docmap){
        $docmap_reverse{$docmap{$_}}=$_;
}


foreach $i (sort {$a <=> $b} keys %docmap){
	print fl_doc "$i:$docmap{$i}\n";
}

close(fl_doc);
open(fl,@ARGV[0]);
@citations=<fl>;
close(fl);
### Produce word mapping
foreach(@citations){
	chomp($_);
        @flds=split(/:/,$_);
        $line="";
  	for($i=2;$i<scalar(@flds);$i++){
    	  $line=join(" ",$line,$flds[$i]);
        }
	@wrds=split(/ /,$line);
	foreach(@wrds){
		$wrd{$_}++;
	}
}

#############################


#############################
foreach $key (sort {$a cmp $b} keys %wrd){
	if($wrd{$key}>40){
		push(@words,$key);
	}
}
open(fl,">wordmap");
$j=0;
%wrdmap=map{$j++ => $_}(@words);
foreach(keys %wrdmap){
	$wrdmap_reverse{$wrdmap{$_}}=$_;
}
foreach $i (sort {$a <=> $b} keys %wrdmap){
	print fl "$i:$wrdmap{$i}:$wrd{$wrdmap{$i}}\n";
}
close(fl);
for($iter=0;$iter<1;$iter++){
#$output="data.lda.train.$iter";
#$input1="citation.train.$iter";
#$input2="content.train.$iter";
#open(fl,">$output");
open(fl,">data");
print fl scalar(keys %docmap)."\n".scalar(keys %wrdmap)."\n".scalar(keys %cited_docs)."\n";
open(citation,"@ARGV[0]");
@citations=<citation>;
close(citation);
foreach(@citations){
	chomp($_);
        @flds=split(/:/,$_);
	if(($flds[0]=~m/[0-9]/)&&($flds[1]=~m/[0-9]/)){
        $line="";
        for($i=2;$i<scalar(@flds);$i++){
          $line=join(" ",$line,$flds[$i]);
        }
        @wrds=split(/ /,$line);
	%temp=();
	foreach (@wrds){
		$temp{$_}++;
	}
	
	foreach (keys %temp){
		if(exists $wrdmap_reverse{$_}){
			print fl "$docmap_reverse{$flds[1]} $docmap_reverse{$flds[0]} $wrdmap_reverse{$_} $temp{$_}\n";
		}
	}
	}
}
open(content,@ARGV[1]);
@content=<content>;
close(content);
foreach(@content){
	chomp($_);
        @flds=split(/:/,$_);
	if($flds[0]=~m/[0-9]/){
        $line="";
        for($i=1;$i<scalar(@flds);$i++){
          $line=join(" ",$line,$flds[$i]);
        }
        @wrds=split(/ /,$line);

	%temp=();
	foreach (@wrds){
		$temp{$_}++;
	}
	#if(not exists $docmap_reverse{$flds[0]}){print "here\n"}
	foreach (keys %temp){
		if(exists $wrdmap_reverse{$_}){
			print fl "$docmap_reverse{$flds[0]} -1 $wrdmap_reverse{$_} $temp{$_}\n";
		}
	}
	}
}
close(fl);
}

for($iter=0;$iter<0;$iter++){
$output="data.lda.test.$iter";
$input1="citation.test.$iter";
$input2="content.test.$iter";
open(fl,">$output");
print fl scalar(keys %docmap)."\n".scalar(keys %wrdmap)."\n1\n";
open(citation,$input1);
@citations=<citation>;
close(citation);
foreach(@citations){
        chomp($_);
        @flds=split(/::/,$_);
        $line="";
        for($i=2;$i<scalar(@flds);$i++){
          $line=join(" ",$line,$flds[$i]);
        }
        @wrds=split(/ /,$line);
        %temp=();
        #$stemmed_words = $stemmer::stem(@wrds);
        #@ids=split(/::/,$citation);
        foreach (@wrds){
                $temp{$_}++;
        }

        foreach (keys %temp){
                if(exists $wrdmap_reverse{$_}){
                        print fl "$docmap_reverse{$flds[1]} $docmap_reverse{$flds[0]} $wrdmap_reverse{$_} $temp{$_}\n";
                }
        }
}        
open(content,$input2);
@content=<content>;
close(content);
foreach(@content){
        chomp($_);
        @flds=split(/::/,$_);
        $line="";
        for($i=1;$i<scalar(@flds);$i++){
          $line=join(" ",$line,$flds[$i]);
        }
        @wrds=split(/ /,$line);

        %temp=();
        foreach (@wrds){
                $temp{$_}++;
        }
        foreach (keys %temp){
                if(exists $wrdmap_reverse{$_}){
                        print fl "$docmap_reverse{$flds[0]} -1 $wrdmap_reverse{$_} $temp{$_}\n";
                }
        }
}
close(fl);
}
