rule centromere_annotation:
    input: assembly="07.Supervised_Final_Assembly/"+config["prefix"]+".assembly.final.fa"
    output: "08.Centromere_Annotation/"+config["prefix"]+".nuclear_genome.centromere.gff3"
    conda: "envs/bioperl_exonerate.yaml"


    shell:
        r"""
        
        
#######################################
# set project-specific variables
        
LRSDAY_HOME=$(realpath LRSDAY)
prefix={config[prefix]} # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.
genome_assembly=$(realpath {input.assembly}) # The path of the input genome assembly.
chrMT_tag={config[chrMT]} # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".
query="$LRSDAY_HOME/data/S288C.centromere.fa" # The S. cerevisiae S288C reference centromere sequences based on Yue et al. (2017) Nature Genetics.
debug={config[debug]} # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".

######################################
# process the pipeline
cd 08.Centromere_Annotation
        
echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m reverse -o $prefix.assembly.nuclear_genome.fa
perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.assembly.nuclear_genome.fa -o $prefix.assembly.nuclear_genome.tidy.fa

exonerate --showvulgar no --showcigar no --showalignment no --showtargetgff yes --bestn 1 $query $prefix.assembly.nuclear_genome.tidy.fa  >$prefix.centromere.exonerate.gff
perl $LRSDAY_HOME/scripts/exonerate_gff2gff3.pl  -i $prefix.centromere.exonerate.gff -o $prefix.centromere.gff3.tmp -t $prefix
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl -r $prefix.assembly.nuclear_genome.tidy.fa -i $prefix.centromere.gff3.tmp -o ../{output} -t $prefix

if [[ $debug == "no" ]]
then
    rm $prefix.assembly.nuclear_genome.fa
    rm $prefix.assembly.nuclear_genome.tidy.fa
    rm $prefix.centromere.exonerate.gff
    rm $prefix.centromere.gff3.tmp
fi

        """
        
