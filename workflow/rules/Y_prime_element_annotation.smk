rule Y_element_annotation:
    conda: "../envs/blat.yaml"
    input: assembly="07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa"
    output:
        gff3="13.Y_Prime_Element_Annotation/"+config["prefix"]+".nuclear_genome.Y_prime_element.gff3",
        fa="13.Y_Prime_Element_Annotation/"+config["prefix"]+".nuclear_genome.Y_prime_element.fa",
        alnfa="13.Y_Prime_Element_Annotation/"+config["prefix"]+".nuclear_genome.Y_prime_element.aln.fa"

    threads: 1

    shell:
        r"""
#######################################
# set project-specific variables
#        
export LRSDAY_HOME=$(realpath LRSDAY)
prefix={config[prefix]} # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.
genome_assembly=$(realpath {input.assembly}) # The path of the input genome assembly.
chrMT_tag={config[chrMT]} # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".
debug={config[debug]}
threads={threads}        
        

        
#######################################
# process the pipeline

## define a safe grep that will not result in pipefail if pattern not found        
c1grep() {{ grep "$@" || test $? = 1; }}

        
cd 13.Y_Prime_Element_Annotation/     

#######################################
# process the pipeline
feature_type="Y_prime_element"
query="$LRSDAY_HOME/data/query.Y_prime_element.long.fa"
length_cutoff_for_completeness=3500

echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m reverse -o $prefix.assembly.nuclear_genome.fa
perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.assembly.nuclear_genome.fa -o $prefix.assembly.nuclear_genome.tidy.fa

blat -maxIntron=1000 $prefix.assembly.nuclear_genome.tidy.fa $query  $prefix.$feature_type.blat.psl
pslCDnaFilter  -minId=0.9 -minAlnSize=1000 -bestOverlap -filterWeirdOverlapped   $prefix.$feature_type.blat.psl  $prefix.$feature_type.blat.filtered.psl
pslScore  $prefix.$feature_type.blat.filtered.psl | sort -nk5 -r > $prefix.$feature_type.blat.filtered.pslScore.out
perl $LRSDAY_HOME/scripts/psl2gff3.pl -i $prefix.$feature_type.blat.filtered.psl -o $prefix.$feature_type.raw.gff3 -ft $feature_type -t $prefix -r $prefix.assembly.nuclear_genome.tidy.fa  -l $length_cutoff_for_completeness
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -i $prefix.$feature_type.raw.gff3 -r $prefix.assembly.nuclear_genome.tidy.fa -o $prefix.nuclear_genome.$feature_type.gff3 -t $prefix
perl $LRSDAY_HOME/scripts/gff2seq_simple.pl -r $prefix.assembly.nuclear_genome.tidy.fa -g $prefix.nuclear_genome.$feature_type.gff3 -o $prefix.nuclear_genome.$feature_type.fa
muscle -in $prefix.nuclear_genome.$feature_type.fa -out $prefix.nuclear_genome.$feature_type.aln.fa

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $prefix.$feature_type.raw.gff3
    rm $prefix.Y_prime_element.blat.psl
    rm $prefix.Y_prime_element.blat.filtered.psl
    rm $prefix.Y_prime_element.blat.filtered.pslScore.out
    rm $prefix.assembly.nuclear_genome.fa
    rm $prefix.assembly.nuclear_genome.tidy.fa
fi

        """
