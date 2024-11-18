#!/usr/bin/env perl

### Script for filtering, joining and re-annotating Ty elements

use strict;
use warnings;
use Carp qw !croak!;
use Array::Utils qw(:all);
use Bio::Tools::RepeatMasker;
use Data::Dumper;
my $debug = 0;

my $gffio = Bio::Tools::GFF->new(-gff_version => 3);
my @Ty = ();

my $maxDiv = 20;
my $maxInternalDist = 100;
## Maximum distance between dissociated Ty fragments to be  joined
my $minLenFracComplete = 0.9; ## Minimum fraction of canonical length to be classified a s complete
my %cannonicalLengths = ( TY1 => 5200,
			  TY2 => 5200,
			  TY3 => 5200,
			  TY4 => 5200,
			  TY5 => 5200,
			  TSU => 5200
			);  
 





my $parser = Bio::Tools::RepeatMasker->new(-file => $ARGV[0]);
### import and filter on sequence divergence

while( my $r = $parser->next_result_2 ) {
    my ($div) = $r->feature1->get_tag_values("perc_div");
    if (ref $r && $div<= $maxDiv) {
	push @Ty, $r;

	print(join "\t", ($r->seq_id(),
		      "perc_div:",($r->feature1->get_tag_values("perc_div")),
		      $r->start,
		      $r->end,
		      $r->strand,
		      $r->score,
		      $r->hstart,
		      $r->hend,
		      $r->hseq_id, ## the TY type
		      $r->primary_tag(),
		      $r->feature1->get_tag_values("myid"),
		      $r->feature1->get_tag_values("is_overlap"),
		      
		     
		      
	      ), "\n") if ($debug);
    }
}

if (@Ty < 1) {
    ### That is not an error per se
    warn "Zero elements after filtering\n";     
    exit 0;
}
print scalar(@Ty)." elements after filtering\n";# if $debug;


sub checkType {
    my @ar = (@_);
    die "at least one element required" unless @ar >0;
    #print @ar,"\n";
    #print @Ty[ar];
    ((grep {ref $_ && $_->hseq_id =~ /T.*-I$/i} @{Ty[@ar]}) > 0) ? "internal" :
	((grep {ref $_ && $_->hseq_id =~ /LTR$/} @{Ty[@ar]}) > 0) ? "soloLTR" : "unknown"
}

sub getTyClass {
    ### internal type takes precedence
    ### If there are different types in one cluster
    ### the returned array has length > 1;
    my @c = map {$_->hseq_id =~ /^(T.*)-I$/i; $1} @{Ty[@_]};
    return unique(@c) if @c > 0;
    @c = map {$_->hseq_id =~ /^(T.*)-LTR$/i; $1} @{Ty[@_]};
    return unique(@c) if @c > 0;
}

sub isComplete {
    return 0 unless @_ > 2;
    my $c = getTyClass(@_)[0]; 
    my $l = (shift) -> start;
    my $r = (pop) -> end;    
    return (($r - $l) >= (cannonicalLength{$c} * $minFractionComplete));
}
    
sub finalSoloType {
    ### simply take the first type if multiple elements
    ### we could check if the LTR types are consistent, but RM does not join
    ### different LTRs
    my $r = $Ty[$_[0]];
    croak $_[0]." is not a valid entry" unless ref $r; 
    $r = $r->hseq_id;
    $r =~ /^T(.+)-LTR$/i;
    die "unknown type $r" unless $1;
    ### Ty1 and Ty2 LTRs are hard to distinguish
    ($1 =~ /^y[12]$/i) ? "TY1/2-soloLTR" : "T$1-soloLTR";
}

sub extendLTR {
    ## we only need to look at left and right element
    ## There is no sanity check if the other internal element(s)
    ## are properly arranged at the moment
    my $l = $_[0]; 
    my $r = $_[@_-1]; ## The same if only one element
    ## we only need to look one element to the left and right
    if ($l > 0 && checkType(($l)) ne 'soloLTR' ) {
	### check if compatible LTR element to the left
	unshift @_, $l-1 if ( ref $Ty[$l-1] &&	    	   
	    $Ty[$l]->seq_id eq $Ty[$l-1]->seq_id &&
	    $Ty[$l]->strand == $Ty[$l-1]->strand &&
	    checkType(($l-1)) eq 'soloLTR' &&
	    ($Ty[$l]->start - $Ty[$l-1]->end < $maxInternalDist)
	    );
    } 
    if ($r < @Ty-1 && checkType(($r)) ne 'soloLTR' ) {
	push @_, $r+1 if (	    	   
	    $Ty[$r]->seq_id eq $Ty[$r+1]->seq_id &&
	    $Ty[$r]->strand == $Ty[$r+1]->strand &&
	    checkType(($r+1)) eq 'soloLTR' &&
	    ($Ty[$r+1]->start - $Ty[$r]->end < $maxInternalDist)
	    );      
    }	
	return @_;
}

## attempt to split custers with more than 1 internal element
## on any internal LTR if present
sub trySplit {
    my @int = grep {$Ty[$_]->} @_;
    return [@_] if @int < 2; # There is only one internal element,
    # no need to split
    # we could split as long as there are at least two internal elements
    while (@int >= 2 &  my $i = shift (@int) ) {
	if ($i < @Ty && checkType(($i+1)) eq 'soloLTR')
	
    }
    
}


### process clusters and join internal or overlapping elements
my @ALL= ();
my $complete = 0;
my $pid = $Ty[0]->feature1->get_tag_values('myid');
# print $Ty[1]->feature1->primary_tag;
my @cluster = (0); 
for (my $i = 1; $i < scalar(@Ty); $i++) {
    #print $Ty[$i]->feature1->get_tag_values('myid');
    print ("processing entry: ",$i,"\n");
  
    my ($cid) = $Ty[$i]->feature1->get_tag_values('myid');
  #  print Dumper $cid;
     #print Dumper "----",$cid,"-----\n";
    # print "currrent cluster: $cid, previous cluster: $pid\n";
     #print Dumper ($Ty[$i]->feature1);
     if ( $cid eq $pid) {	
	 push  @cluster, $i;	
	 next;
     } else {
	 ## Time to start a new cluster:
	 ## And to Analyze the current one    
     
	 ### The consecutive cluster ids are not the same 
	 #next ELEM if (@cluster < 1);
	 ### Check if the cluster is complete
	 ### check if it contains an Internal element
	 ### Clusters without internal element are never complete and are
	 ### re-annotated as SoloLTR, however many there are
	 die "empty cluster encountered, this should never happen!"
	     unless @cluster > 0;    	 
	 my $type = checkType(@cluster);
	 ### less than 3 elements are never complete
	 ### But does it contain internal elements?
	 ### cluster with 3 elements may be correct but must be of
	 ### correct syntax TyX-LTR -> TyX-I -> TyX-LTR
	 ### Because extendLTR checks the grammar, we can use it
	 ### all clusters using the same logic
	 if ($type eq "soloLTR") {
	     print @cluster,"\n";
	     $type = finalSoloType(@cluster);
	 } elsif ($type eq "internal") {
	     ### try to extend the cluster with LTR on both sides if missing
	     @cluster = extendLTR(@cluster);
	 }
	 if (scalar(@cluster < 3)){		 
	  
	     print "cluster $pid of type $type  with length ".scalar(@cluster)."\n" ;	     
	 }
	 elsif (@cluster == 3) {

	 print "Potentially Complete cluster$pid  of type $type  with length == 3\n";
	 ### At this point we just need to check whether we have to extend
     } else {
	 ### That may be an overlapping cluster or contain multiple HSPs
	 print "Potentially complete cluster $pid of type $type with length > 3, trying to split\n";
	 ###
	 my @clusters = trySplit(@cluster);


    }

     
	 ($pid) = $Ty[$i]->feature1()->get_tag_values('myid');
	 push @ALL, [@cluster]; # save the current cluster for the final processing step
	 @cluster = ($i);


     }

}












package Bio::Tools::RepeatMasker;
no warnings qw(redefine);

sub next_result_2 {
    my ($self) = @_;
    
    local $_;
    while (defined($_=$self->_readline()) ) {
	if (/no repetitive sequences detected/) {
	    $self->warn( "RepeatMasker didn't find any repetitive sequences\n");
	    return ;
	}
	#ignore introductory lines
	if (/\d+/) {
	    my @element = split;
	    # ignore features with negatives
	    next if ($element[11-13] =~ /-/);
	    my (%feat1, %feat2);
	    my @line = split;
	    my ($score, $query_name, $query_start, $query_end, $strand,
		$repeat_name, $repeat_class ) = @line[0, 4, 5, 6, 8, 9, 10];
	    my ($perc_div, $perc_del, $perc_ins, $id, $is_overlap) = @line[1..3,14,15];

	    my ($hit_start,$hit_end);

	    if ($strand eq '+') {
		($hit_start, $hit_end) = @line[11, 12];
		$strand = 1;
	    } elsif ($strand eq 'C') {
		($hit_end, $hit_start) = @line[12, 13];
		$strand = -1;
	    }
	    my $rf = Bio::SeqFeature::Generic->new
		(-seq_id      => $query_name,
		 -score       => $score,
		 -start       => $query_start,
		 -end         => $query_end,
		 -strand      => $strand,
		 -source_tag  => 'RepeatMasker',
		 -primary_tag => $repeat_class,
		 -tag => {
		     'Target'=> [$repeat_name, $hit_start, $hit_end],
			 'perc_div' => $perc_div,
			 'perc_del' => $perc_del,
			 'perc_ins' => $perc_ins,
			 'myid' => [$id],
			 'is_overlap' =>
			 ($is_overlap && $is_overlap eq "*") ? 1 : 0 

		 },
		);

	    my $rf2 = Bio::SeqFeature::Generic->new
		(-seq_id         => $repeat_name,
		 -score          => $score,
		 -start          => $hit_start,
		 -end            => $hit_end,
		 -strand         => $strand,
		 -source_tag     => "RepeatMasker",
		 -primary_tag    => $repeat_class,
		 -tag => { 'Target'=> [$query_name,$query_start,$query_end] },
		);

	    my $fp = Bio::SeqFeature::FeaturePair->new(-feature1 => $rf,
						       -feature2 => $rf2);
	    return $fp;
	}
    }
}
