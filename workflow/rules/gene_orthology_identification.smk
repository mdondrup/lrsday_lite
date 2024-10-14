prefix=config["prefix"]

rule gene_orthology_identification:
    conda: "envs/proteinortho.yaml"
    input: 
           input_nuclear_genome_gff="09.Nuclear_Gene_Annotation/"+prefix+".nuclear_genome.gff3",
           query_nuclear_genome_PoFF_faa="09.Nuclear_Gene_Annotation/"+prefix+".nuclear_genome.PoFF.faa", 
           query_nuclear_genome_PoFF_gff="09.Nuclear_Gene_Annotation/"+prefix+".nuclear_genome.PoFF.gff", 

           input_mitochondrial_genome_gff="10.Mitochondrial_Gene_Annotation/"+prefix+".mitochondrial_genome.gff3",
           query_mitochondrial_genome_PoFF_faa="10.Mitochondrial_Gene_Annotation/"+prefix+".mitochondrial_genome.PoFF.faa",
           query_mitochondrial_genome_PoFF_gff="10.Mitochondrial_Gene_Annotation/"+prefix+".mitochondrial_genome.PoFF.gff",

    output:
        nuc_tsv="14.Gene_Orthology_Identification/"+prefix+".nuclear_genome.proteinortho.tsv",
        mit_tsv="14.Gene_Orthology_Identification/"+prefix+".mitochondrial_genome.proteinortho.tsv"
    log:
        nuc="14.Gene_Orthology_Identification/"+prefix+".nuclear_genome.proteinortho.log",
        mit="14.Gene_Orthology_Identification/"+prefix+".mitochondrial_genome.proteinortho.log"

    threads: 80

    shell:
        r"""
#######################################
# set project-specific variables
#        
export LRSDAY_HOME=$(realpath LRSDAY)
prefix={config[prefix]} # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.

chrMT_tag={config[chrMT]} # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".
debug={config[debug]}
threads={threads}        
        

        
#######################################

ref_PoFF_faa="$LRSDAY_HOME/data/SGDref.PoFF.faa" # The file path of the reference proteome file in FASTA format: for S. cerevisiae and its close relatives, you can directly use the pre-shipped file: SGDref.PoFF.faa; if you work with other organisms, you can check ProteinOrtho's manual for details on how to prepare such file.
ref_PoFF_gff="$LRSDAY_HOME/data/SGDref.PoFF.gff" # The path of the reference gene GFF file in GFF format: for S. cerevisiae and its close relatives, you can directly use the pre-shipped file: SGDref.PoFF.gff; if you work with other organisms, you can check ProteinOrtho's manual for details on how to prepare such file.



#######################################
# process the pipeline

cd 14.Gene_Orthology_Identification/

cp $ref_PoFF_faa ref.PoFF.faa
cp $ref_PoFF_gff ref.PoFF.gff

## sanitize input files:

sed -i -E '/^>/! s/[^XOUBZACDEFGHIKLMNPQRSTVWYxoubzacdefghiklmnpqrstvwy]//g; /^$/d' ref.PoFF.faa;

        
# process nuclear gene annotation

    
    cp -v ../{input.query_nuclear_genome_PoFF_faa} query.nuclear_genome.PoFF.faa
    cp -v ../{input.query_nuclear_genome_PoFF_gff} query.nuclear_genome.PoFF.gff
    proteinortho_nuclear_genome_input="query.nuclear_genome.PoFF.faa ref.PoFF.faa"
    proteinortho -cpus=$threads -p=blastp -binpath=$CONDA_PREFIX/bin -singles -synteny -project="$prefix.nuclear_genome" $proteinortho_nuclear_genome_input > ../{log.nuc} 2>&1
    perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i ../{input.input_nuclear_genome_gff} -x $prefix.nuclear_genome.proteinortho.tsv -r ref.PoFF.faa -q query.nuclear_genome.PoFF.faa -o $prefix.nuclear_genome.SGD_orthology_mapped.gff3


    if [[ $debug == "no" ]]
    then
	rm $prefix.nuclear_genome.ffadj-graph
	rm $prefix.nuclear_genome.blast-graph
	rm $prefix.nuclear_genome.proteinortho-graph
	rm $prefix.nuclear_genome.poff-graph
    fi


# process mitochondrial gene annotation

   
    cp -v ../{input.query_mitochondrial_genome_PoFF_faa} query.mitochondrial_genome.PoFF.faa
    cp -v ../{input.query_mitochondrial_genome_PoFF_gff} query.mitochondrial_genome.PoFF.gff
    proteinortho_mitochondrial_genome_input="query.mitochondrial_genome.PoFF.faa ref.PoFF.faa"
    proteinortho -cpus=$threads -p=blastp -binpath=$CONDA_PREFIX/bin -singles -synteny -project="$prefix.mitochondrial_genome"  $proteinortho_mitochondrial_genome_input > ../{log.mit} 2>&1
    perl $LRSDAY_HOME/scripts/update_gff3_by_proteinortho.pl -i ../{input.input_mitochondrial_genome_gff} -x $prefix.mitochondrial_genome.proteinortho.tsv -r ref.PoFF.faa -q query.mitochondrial_genome.PoFF.faa -o $prefix.mitochondrial_genome.SGD_orthology_mapped.gff3

    if [[ $debug == "no" ]]
    then
	rm $prefix.mitochondrial_genome.ffadj-graph
	rm $prefix.mitochondrial_genome.blast-graph
	rm $prefix.mitochondrial_genome.proteinortho-graph
	rm $prefix.mitochondrial_genome.poff-graph
    fi



# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm query.*
    rm ref.*
    rm proteinortho_cache_${{prefix}}.log
fi
    





        

        """
