The script is designed to process the output of a RepeatMasker run, focusing on filtering, clustering, joining, and re-annotating Ty elements (a type of retrotransposon commonly studied in yeast genomes). Below is an explanation of its key functions:

= Core Functionality =

Input Parsing

The script reads RepeatMasker results from a file provided as a command-line argument.
It uses the Bio::Tools::RepeatMasker module to parse and extract information about each repeat feature, including:
Sequence ID
Divergence percentage
Start and end positions
Repeat type (e.g., LTR, internal)
Filtering

Filters repeat elements based on divergence percentage (perc_div) using a maximum divergence threshold ($maxDiv).
Only repeat elements with a divergence less than or equal to this threshold are processed further.
Clustering

Groups Ty-related elements into clusters based on shared identifiers (myid) and proximity in the genome.
Handles overlapping or nearby repeat elements and attempts to join fragments into cohesive clusters.
Classification

Each cluster is classified into one of the following types:
Internal: Internal sequence of the retrotransposon.
Solo LTR: Long Terminal Repeat (LTR) without internal components.
Complete Ty Element: Includes LTRs flanking internal regions, meeting the criteria for completeness.
Cluster Analysis and Annotation

Determines if clusters represent:
Complete elements: Based on predefined canonical lengths and configuration.
Truncated elements: If the cluster doesn't meet completeness criteria.
Updates annotations for each cluster to indicate:
Type (e.g., TY1-LTR, TY1/2-soloLTR, TY3-I)
Completeness (yes/no)
Cluster ID
Extends clusters by incorporating adjacent LTRs when appropriate.
Re-Annotation

Modifies repeat feature annotations to reflect:
Cluster type
Completeness status
Updated names (e.g., TySon.TY1_1, TySon.TY1/2-soloLTR_2).
Output

Outputs:
BED file: Genome coordinates of the processed repeat clusters.
Counts file: Statistics on the number of elements of each type, truncated and complete.
(Optionally) GFF file: Re-annotated features in GFF format.
Key Features
Handles Fragmented Ty Elements

Merges disconnected or overlapping fragments to recover complete retrotransposons.
Customizable Parameters

maxDiv: Maximum allowed sequence divergence for filtering.
maxInternalDist: Maximum allowed distance between fragments to join them.
minFractionComplete: Minimum fraction of canonical length to qualify as complete.
Support for Ty Families

Includes canonical lengths for different Ty families (e.g., TY1, TY2, TY3).
Adjusts annotations accordingly when multiple types are present in a cluster.
Debugging and Output Options

Debug mode provides detailed logs.
Flexible output formats: BED, GFF, and counts.
Intended Usage
This script is typically used in genome annotation pipelines to refine and classify retrotransposons after an initial RepeatMasker run.
Particularly suited for studying yeast Ty elements, but adaptable to similar retrotransposons in other organisms.
