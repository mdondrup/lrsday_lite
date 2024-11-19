#!/usr/bin/env perl


## Script for filtering, joining, and re-annotating Ty elements

use strict;
use warnings;
use Carp qw(croak carp confess);
use Getopt::Long;
use Pod::Usage;
use Array::Utils qw(:all);
use Bio::Tools::RepeatMasker;
use Data::Dumper;

### Default parameters
my $debug = 1;
my $bed = 1;
my $counts = 1;
my $prefix = 'TySon';
my $class_filter = '';
my $complete_only = 0;
my $maxDiv = 30;            # Maximum divergence percentage
my $maxInternalDist = 100;  # Maximum distance to join fragments
my $minFractionComplete = 0.9; # Minimum fraction of canonical length

### Add command-line options
my $help = 0;    # Help flag
GetOptions(
    'maxDiv=f'             => \$maxDiv,              # Max divergence (float)
    'maxInternalDist=i'    => \$maxInternalDist,      # Max internal distance (integer)
    'minFractionComplete=f'=> \$minFractionComplete,  # Min fraction complete (float)
    'debug=i'        => \$debug,         # Debug level
    'bed!'           => \$bed,          # Enable/disable BED output
    'counts!'        => \$counts,       # Enable/disable counts output
    'prefix=s'       => \$prefix,       # Prefix for names
    'class_filter=s' => \$class_filter, # Filter for specific class
    'complete_only!' => \$complete_only, # Output only complete clusters
    'help|h'         => \$help,         # Help flag
) or pod2usage(2);

pod2usage(1) if $help;




###############################################################
################## Global variables ###########################

my $gffio = Bio::Tools::GFF->new(-gff_version => 3); ## gff version to use
my @Ty = ();
my %counts = ();

my $maxDiv = 20;
my $maxInternalDist = 100;
## Maximum distance between dissociated Ty fragments to be  joined
my $minFractionComplete = 0.9; ## Minimum fraction of canonical length to be classified a s complete

my %cannonicalLengths = ( 'TY1' => 5924,
			  'TY2' => 5959,
			  'TY3' => 5351,
			  'TY3_1p' => 5360,
			  'TY4' => 6223,
			  'TY5' => 5349,
			  'TSU4' => 5997
			);  


my $parser = Bio::Tools::RepeatMasker->new(-file => $ARGV[0]);
### import and filter on sequence divergence

while( my $r = $parser->next_result_2 ) {
    my ($div) = $r->feature1->get_tag_values("perc_div");
    if (ref $r && $div<= $maxDiv) {
	push @Ty, $r;
## print the parsed input
	print STDERR (join "\t", ($r->seq_id(),
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
	      ), "\n") if ($debug > 1);
    }
}

### END imported the input

if (@Ty < 1) {
    ### That is not an error per se
    carp "Zero elements after filtering for percent divergence\n";     
    exit 0;
}
print STDERR  scalar(@Ty)." elements after filtering\n" if $debug;


sub checkType {
     confess "call to checkType contained undef!" if ((grep { (! defined ($_)) } @_) > 0);
    confess "at least one element required for checkType " unless @_ > 0;  
    
    ((grep {ref $_ && $_->hseq_id =~ /T.*-I$/i} @{Ty[@_]}) > 0) ? "internal" :
	((grep {ref $_ && $_->hseq_id =~ /LTR$/} @{Ty[@_]}) > 0) ? "soloLTR" : "unknown"
}

sub getTyClass {
    ### internal type takes precedence
    ### If there are different types in one cluster
    ### the returned array has length > 1;
    my @c = map {$_->hseq_id =~ /^(T.*)-I$/i; uc($1) if $1} @{Ty[@_]};
    return (unique(@c, @c)) if @c > 0;
    @c = map {$_->hseq_id =~ /^(T.*)-LTR$/i;  uc($1) if $1} @{Ty[@_]};
    return (unique(@c, @c)) if @c > 0;
    #confess "unknown Ty class";  
    return ();
}

sub getScoreSum {
    # calculate the combined score for the cluster
    my $score = 0;
    map {$score += $_->score} @{Ty[@_]};
    return $score;
}


sub isComplete {
    
    @_ = grep {defined $_} @_; 
    return 0 unless @_ > 2;
    ## check if both ends are LTRs;
    return 0 unless (checkType($_[0]) eq 'soloLTR' && checkType($_[-1]) eq 'soloLTR' );
    my @ar = grep {$_} getTyClass(@_);
    my $c = join '/', @ar;
    #print STDERR "----$c";
    my $l = @Ty[(shift)] -> start;
    my $r = @Ty[(pop)] -> end;
    confess "unknown Ty type $c" unless exists $cannonicalLengths{$c};   
    return (($r - $l) >= $cannonicalLengths{$c} * $minFractionComplete) ? 1 : 0;
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
    my $r = $_[-1]; ## The same if only one element
    confess "$r isn't defined " unless defined $r;
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
    if ( $r < @Ty-1 && checkType(($r)) ne 'soloLTR' ) {
	push @_, $r+1 if (	    	   
	    $Ty[$r]->seq_id eq $Ty[$r+1]->seq_id &&
	    $Ty[$r]->strand == $Ty[$r+1]->strand &&
	    checkType(($r+1)) eq 'soloLTR' &&
	    ($Ty[$r+1]->start - $Ty[$r]->end < $maxInternalDist)
	    );      
    }	
	return (@_);
}

## attempt to split custers with more than 1 internal element
## on any internal LTR if present
sub trySplit {
    my @int = grep { checkType($_) eq 'internal' } @_;
    @_ = grep { defined $_ } @_;
    return [@_] if @int < 2; # There is only one internal element,
    # no need to split
    # we could split as long as there are at least two internal elements
    my $last = 0; ## store the first element for a new cluster, initially that's 0
    my @ret =();
    foreach my $i (0..scalar(@_-1)) {
	last unless @int >=2; # no point splitting further if only 1 left
	shift (@int) if checkType(@_[$i]) eq 'internal';
	if ($i+1 < @_ && checkType(@_[$i]) eq 'internal' &&
	    checkType(@_[$i+1]) eq 'soloLTR') {
	    push @ret, [@_[$last..$i+1]];	   
	    if ($i+2 < @_ && checkType($i+2) eq 'soloLTR') {
		$last = $i+2;
	    } else {
		$last = $i+1; # There's only a single LTR betwen two I's
	    }	    
	}
    }
    # add the rest of the cluster, if any. this includes the case of no hit at all
    if ($last < @_) {
	push @ret, [@_[$last..@_-1]];
    }
    return (@ret);
}


### process clusters and join internal or overlapping elements
my @ALL= ();
my $pid = $Ty[0]->feature1->get_tag_values('myid');
my @cluster = (0); 
foreach my $i (1..scalar(@Ty)-1) {
   
    print STDERR  ("processing entry: ",$i,"\n") if $debug;
  
    my ($cid) = $Ty[$i]->feature1->get_tag_values('myid');
 
    if ( $cid eq $pid) {
	 push  (@cluster, $i);	
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
	 my $type = checkType((@cluster));
	 if ($debug) {
	     print STDERR  "Cluster syntax before editing:\t\t";
	     print STDERR  (join ',', (map {$Ty[$_]->hseq_id} @cluster)); print STDERR  "\n" 
	 }
	 ### less than 3 elements are never complete
	 ### But does it contain internal elements?
	 ### cluster with 3 elements may be correct but must be of
	 ### correct syntax TyX-LTR -> TyX-I -> TyX-LTR
	 ### Because extendLTR checks the grammar, we can use it
	 ### all clusters using the same logic
	 if ($type eq "soloLTR") {	    
	     $type = finalSoloType(@cluster);
	 } elsif ($type eq "internal") {
	     ### try to extend the cluster with LTR on both sides if missing
	     @cluster = extendLTR(@cluster);
	     if ($debug) {
		 print STDERR  "Cluster syntax after LTR extension:\t";
		 print STDERR  (join ',', (map {$Ty[$_]->hseq_id} @cluster)); print STDERR  "\n" 
	     }
	 }
	 if (scalar(@cluster < 3)){		 
	  
	     print STDERR  "cluster $pid of type $type  with length ".scalar(@cluster)."\n" if $debug ;
	     push @ALL, [@cluster]; # save the current cluster for the final processing step
	 }
	 elsif (@cluster == 3) {
	     print STDERR  "Potentially Complete cluster$pid  of type $type  with length == 3\n" if $debug;
	     push @ALL, [@cluster]; # save the current cluster for the final processing step
     } else {
	 ### That may be an overlapping cluster or contain multiple HSPs
	 print STDERR  "Potentially complete cluster $pid of type $type with length > 3, trying to split\n" if $debug;
	 ###
	 @cluster = trySplit(@cluster);
	 if (@cluster >1 && $debug) {
	     print STDERR  "--> cluster was split into ". scalar(@cluster)," clusters\n";
	     print STDERR  join ' ', ( map { '[' . (join ',', (map {$Ty[$_]->hseq_id} @$_)) . ']' } @cluster );
	 }
	 push @ALL, @cluster;
    }
	 ($pid) = $Ty[$i]->feature1()->get_tag_values('myid');	
	 @cluster = ($i);
     }    
} ### END processing of clusters

#### print all cluster elements in the way they are after modification

foreach my $c (0..@ALL-1) {
   
    my @cluster = @{$ALL[$c]};
    my $class = '';
    my $complete = 0;
    my $name = $prefix.'.';
    if ( checkType(@cluster) eq 'soloLTR' ) {
	$class = finalSoloType(@cluster);
	next if $complete_only;
	next if  ($class_filter && $class !~ /^$class_filter/);
	map {$Ty[$_]->feature1->add_tag_value(Class=>$class);
	     $Ty[$_]->feature1->add_tag_value(complete=>'No')} @cluster;
	$name .= $class."_".($c+1);
	$counts{soloLTR}++;
	$counts{$class}++;
    } elsif ( checkType(@cluster) eq 'internal' ) {
	$class = join( '/', (grep {$_} (getTyClass(@cluster)))); ## accomodate for divergent classes, not so great but rarely happens
        $complete = isComplete(@cluster);
	map {$Ty[$_]->feature1->add_tag_value(Class=> $class . ( (checkType($_) eq 'soloLTR' )? '-LTR' : '-I' ))} @cluster;
	map {$Ty[$_]->feature1->add_tag_value(complete=> ($complete)?'Yes':'No')} @cluster;
	$name .=  $class."_".($c+1).(($complete) ? '' : '_truncated');
	$counts{total_truncated}++ && 	$counts{$class."_truncated"}++ if not $complete;
	$counts{total_complete}++  && 	$counts{$class."_complete"}++  if $complete;
    } else {
	confess "unknown cluster type: " .  checkType(@cluster);
    }
    next if ($class_filter && $class ne $class_filter);
    next if ($complete_only && ! $complete);
   
    
    
    
    #print STDERR $c;
    map {$Ty[$_]->feature1->add_tag_value(cluster => $c+1)} @cluster;
    map {$Ty[$_]->feature1->add_tag_value(Name => $name)} @cluster;
    if ($bed) {
	### write a bed file with genome coordinates for the full matches
	my $l = $Ty[$cluster[0]];
	my $r = $Ty[$cluster[-1]];
	print join "\t", ($l->seq_id(),
			  $l->start,
			  $r->end,
			   $name,
			  getScoreSum(@cluster),
			  ($r->strand > 0)? '+' : '-'
	    );
	print "\n";		     
	
    } else {
	map {print $Ty[$_]->feature1->gff_string($gffio); print "\n"} @cluster;
    }
}

if ($counts) {
    *COUNTS = *STDOUT; # default is to also print ot STDOUT if no prefix specified
    open COUNTS, ">$prefix.counts" || die "couldn't open counts file: $!" if ($prefix); 
    my @k = sort(keys %counts);
    print COUNTS (join "\t", @k ),"\n";
    print COUNTS (join "\t", @counts{@k}), "\n";
    #print Dumper [\%counts];
}



##########################################################################################################

package Bio::Tools::RepeatMasker;
no warnings qw(redefine);

### The original function doesn't read the cluster ids assigned by RM

sub next_result_2 {
    my ($self) = @_;
    
    local $_;
    while (defined($_=$self->_readline()) ) {
	if (/no repetitive sequences detected/) {
	    $self->warn( "RepeatMasker didn't find any repetitive sequences\n");
	    return ;
	}
	next if /^\s*$/;
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


__END__

=head1 NAME

TySon.pl - Filter, join, and re-annotate Ty elements.

=head1 SYNOPSIS

TyFilter [options] input_file

 Options:   
   --maxDiv                 Set maximum sequence divergence (default: 20)
   --maxInternalDist        Set internal distance (default: 100)
   --minFractionComplete    Set fraction to call a cluster complete (default: 0.9)
   --debug=<level>          Set debug level (default: 1)
   --[no]_bed                Enable BED file output (default: enabled)
   --counts                 Enable counts file output (default: enabled)
   --prefix=<string>        Prefix for names (default: "TySon")
   --class_filter=<regex>   Filter for a specific class of Ty elements
   --complete_only          Output only complete clusters
   --help                   Display this help and exit

=head1 DESCRIPTION

This script processes RepeatMasker output to filter, join, and re-annotate Ty elements, a type of retrotransposon commonly studied in yeast genomes, based on various criteria such as sequence divergence and canonical length.

=head1 OPTIONS

=over 4

=item B<--maxDiv>

Set maximum divergence (float). Filters on RepeatMasker output. [20]

=item B<--maxInternalDist>

Set maximum internal distance for soloLTR elements to be joined with incomplete clusters (integer). [100]

=item B<--minFractionComplete>

Set minimum fraction of canonical Ty-element length to classify as complete(float). Filters on RepeatMasker output. [0.9]

=item B<--debug>

Set the debug level (1,2). Higher values provide more verbose output.

=item B<--bed, --no_bed>

Enable or disable BED file output. Enabled by default. GFF3 Otherwise

=item B<--counts>, B<--no-counts>

Enable or disable counts file output. Enabled by default.

=item B<--prefix>

Specify the prefix to use for cluster names. Default is "TySon".

=item B<--class_filter>

Provide a regex to filter for specific Ty element classes.

=item B<--complete_only>

If specified, outputs only complete clusters. By default, both complete and incomplete clusters are output.

=item B<--help>

Display a brief help message and exit.

=back

=head1 EXAMPLES

Run with default options:

  perl TySon.pl input_file

Run with custom prefix and BED file disabled:

  perl TySon.pl --prefix="CustomPrefix" --maxDiv 25 --maxInternalDist 50 --minFractionComplete 0.8 --debug 1


=head1 AUTHOR

Your Name <michael.dondrup <at> uib.no>

=cut
