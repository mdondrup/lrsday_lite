rule X_element_annotation:
    conda: "envs/mfannot.yaml"
    input: assembly="07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa"
    output:
        gff3="12.Core_X_Element_Annotation/"+config["prefix"]+".nuclear_genome.X_element.gff3",
        fa="12.Core_X_Element_Annotation/"+config["prefix"]+".nuclear_genome.X_element.fa",
        alnfa="12.Core_X_Element_Annotation/"+config["prefix"]+".nuclear_genome.X_element.aln.fa"

    threads: 80

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

        
cd 12.Core_X_Element_Annotation/     

feature_type="X_element"
length_cutoff_for_completeness=300

echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m reverse -o $prefix.assembly.nuclear_genome.fa
perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.assembly.nuclear_genome.fa -o $prefix.assembly.nuclear_genome.tidy.fa

nhmmer -E 1 --tblout $prefix.$feature_type.nhmmer.out $LRSDAY_HOME/data/S288C.$feature_type.hmm $prefix.assembly.nuclear_genome.tidy.fa
perl $LRSDAY_HOME/scripts/nhmer2seq.pl -r $prefix.assembly.nuclear_genome.tidy.fa -i $prefix.$feature_type.nhmmer.out -e 0.0001 -p $prefix -ft $feature_type -l $length_cutoff_for_completeness
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -i $prefix.$feature_type.raw.gff3 -r $prefix.assembly.nuclear_genome.tidy.fa -o $prefix.nuclear_genome.$feature_type.gff3 -t $prefix
perl $LRSDAY_HOME/scripts/gff2seq_simple.pl -r $prefix.assembly.nuclear_genome.tidy.fa -g $prefix.nuclear_genome.$feature_type.gff3 -o $prefix.nuclear_genome.$feature_type.fa 
muscle -align $prefix.nuclear_genome.$feature_type.fa -output $prefix.nuclear_genome.$feature_type.aln.fa

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $prefix.$feature_type.raw.gff3
    rm $prefix.$feature_type.raw.fa
    rm $prefix.assembly.nuclear_genome.fa
    rm $prefix.assembly.nuclear_genome.tidy.fa
fi
        """        
        
