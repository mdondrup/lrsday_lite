prefix=config["prefix"]

rule annotation_integration:
    conda: "envs/bioperl_exonerate.yaml"
    input:
        genome_assembly="07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa",              
        centromere_gff3="08.Centromere_Annotation/"+prefix+".nuclear_genome.centromere.gff3",
        TE_gff3="11.TE_Annotation/"+prefix+".nuclear_genome.TE.gff3", 
        X_element_gff3="12.Core_X_Element_Annotation/"+prefix+".nuclear_genome.X_element.gff3",
        Y_prime_element_gff3="13.Y_Prime_Element_Annotation/"+prefix+".nuclear_genome.Y_prime_element.gff3", 
        nuclear_gene_gff3="14.Gene_Orthology_Identification/"+prefix+".nuclear_genome.SGD_orthology_mapped.gff3", 
        mitochondrial_gene_gff3="14.Gene_Orthology_Identification/"+prefix+".mitochondrial_genome.SGD_orthology_mapped.gff3"
           

    output:
        nuc_final_fa="15.Annotation_Integration/"+prefix+".nuclear_genome.tidy.fa",
        nuc_final_gff3="15.Annotation_Integration/"+prefix+".nuclear_genome.tidy.gff3",
        nuc_final_cds_fa="15.Annotation_Integration/"+prefix+".nuclear_genome.tidy.cds.fa",
        mit_final_fa="15.Annotation_Integration/"+prefix+".mitochondrial_genome.tidy.fa",
        mit_final_gff3="15.Annotation_Integration/"+prefix+".mitochondrial_genome.tidy.gff3",
        mit_final_cds_fa="15.Annotation_Integration/"+prefix+".mitochondrial_genome.tidy.cds.fa"   
    params:
        chrMT_genetic_code_table=3
    threads: 80

    shell:
        r"""
#######################################
# set project-specific variables
#        
export LRSDAY_HOME=$(realpath LRSDAY)
prefix={config[prefix]} # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.
genome_assembly=$(realpath {input.genome_assembly}) # The path of the input genome assembly.
chrMT_tag={config[chrMT]} # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".
debug={config[debug]}
threads={threads}

######################################
# process the pipeline
set -x
cd 15.Annotation_Integration/        
echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m reverse -o $prefix.nuclear_genome.fa
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m normal -o $prefix.mitochondrial_genome.fa

perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.nuclear_genome.fa -o $prefix.nuclear_genome.tidy.fa

cat ../{input.centromere_gff3} ../{input.TE_gff3} ../{input.X_element_gff3} ../{input.Y_prime_element_gff3} ../{input.nuclear_gene_gff3}  >$prefix.nuclear_genome.concatenated.gff3
perl $LRSDAY_HOME/scripts/sort_gff3.pl -i $prefix.nuclear_genome.concatenated.gff3 -t $prefix -o $prefix.nuclear_genome.tidy.gff3 -r $prefix.nuclear_genome.tidy.fa


    perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl -r $prefix.nuclear_genome.tidy.fa -g $prefix.nuclear_genome.tidy.gff3 -o $prefix.nuclear_genome.tidy.cds.fa
    perl $LRSDAY_HOME/scripts/cds2protein.pl -i $prefix.nuclear_genome.tidy.cds.fa -p $prefix.nuclear_genome.tidy -t 1


    perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.mitochondrial_genome.fa -o $prefix.mitochondrial_genome.tidy.fa
    cp ../{input.mitochondrial_gene_gff3} $prefix.mitochondrial_genome.tidy.gff3
  
	perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl -r $prefix.mitochondrial_genome.tidy.fa -g $prefix.mitochondrial_genome.tidy.gff3 -o $prefix.mitochondrial_genome.tidy.cds.fa
	perl $LRSDAY_HOME/scripts/cds2protein.pl -i $prefix.mitochondrial_genome.tidy.cds.fa -p $prefix.mitochondrial_genome.tidy -t {params.chrMT_genetic_code_table}
    



# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $prefix.nuclear_genome.concatenated.gff3

    if [[ -e $prefix.nuclear_genome.fa ]]
    then
	rm $prefix.nuclear_genome.fa
    fi

    if [[ -e $prefix.mitochondrial_genome.fa ]]
    then
	rm $prefix.mitochondrial_genome.fa
    fi
fi

        """
        
