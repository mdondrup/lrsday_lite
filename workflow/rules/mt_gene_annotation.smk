rule mt_gene_annotation:
    conda: "../../envs/mfannot.yaml"
    input:
        assembly="07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa"
    output:
         gff3= "10.Mitochondrial_Gene_Annotation/"+config["prefix"]+".mitochondrial_genome.gff3",
         cdsfa= "10.Mitochondrial_Gene_Annotation/"+config["prefix"]+".mitochondrial_genome.cds.fa",
         poffgff= "10.Mitochondrial_Gene_Annotation/"+config["prefix"]+".mitochondrial_genome.PoFF.gff",
         poffffn= "10.Mitochondrial_Gene_Annotation/"+config["prefix"]+".mitochondrial_genome.PoFF.ffn",
         pofffaa= "10.Mitochondrial_Gene_Annotation/"+config["prefix"]+".mitochondrial_genome.PoFF.faa"         
    shell:
        r"""
#######################################
# load environment variables for LRSDAY

export LRSDAY_HOME=$(realpath LRSDAY)
prefix={config[prefix]} # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.
genome_assembly=$(realpath {input.assembly}) # The path of the input genome assembly.
chrMT_tag={config[chrMT]} # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".
debug={config[debug]} 
        
source tools/bash.env
        
export PERL5LIB="$pirobject_dir/lib"
export RNAFINDER_CFG_PATH="$rnafinder_dir"
export MF2SQN_LIB="$mf2sqn_dir/lib"
export MFANNOT_LIB_PATH="$mfannot_data_dir/protein_collections"
export MFANNOT_EXT_CFG_PATH="$mfannot_data_dir/config"
export MFANNOT_MOD_PATH="$mfannot_data_dir/models"
export BLASTMAT="$blast_matrices_dir"
export EGC="$mfannot_data_dir/EGC"
export ERPIN_MOD_PATH="$mfannot_data_dir/models/Erpin_models"
export PIR_DATAMODEL_PATH="$pirobject_dir/PirModels"
export PATH="$flip_dir:$umac_dir:$erpin_dir/bin:$pirobject_dir:$pirmodels_dir:$hmmsearchwc_dir:$mf2sqn_dir:$mf2sqn_dir:$grab_fasta_dir:$rnafinder_dir:$mfannot_dir:$PATH"

export trnascan_dir=$CONDA_PREFIX/lib/tRNAscan-SE

        
#######################################
# set project-specific variables

genetic_code_table=3 # The NCBI genetic code table (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for the annotated mitochondrial genome. Default = 3 (i.e. Yeast Mitochondria)

######################################
# process the pipeline
set -x
cd 10.Mitochondrial_Gene_Annotation
echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list

perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m normal -o $prefix.assembly.mitochondrial_genome.fa

if [[ $(egrep -c "^>" "$prefix.assembly.mitochondrial_genome.fa") -eq 0 ]]
then
    echo "No mitochondrial genome assembly was detected! Check MT prefix!"
    exit 1    
else
    
    perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.assembly.mitochondrial_genome.fa -o $prefix.assembly.mitochondrial_genome.tidy.fa

    mkdir -p $$.tmp
   ## mfannot uses #!/usr/bin/perl     
   perl $mfannot_dir/mfannot \
	--genetic $genetic_code_table \
	--outputfile $prefix.mitochondrial_genome.mfannot.out.txt \
	--logfile $prefix.mitochondrial_genome.mfannot.log \
	--T $(pwd)/$$.tmp \
	$prefix.assembly.mitochondrial_genome.tidy.fa

    tRNAscan-SE -D -g $trnascan_dir/gcode.ystmito $prefix.assembly.mitochondrial_genome.tidy.fa  -o $prefix.assembly.mitochondrial_genome.trnascan.out.txt


    perl $LRSDAY_HOME/scripts/mfannot2gff3.pl -mfannot_out $prefix.mitochondrial_genome.mfannot.out.txt -trnascan_out $prefix.assembly.mitochondrial_genome.trnascan.out.txt -o $prefix.mitochondrial_genome.pre_sort.gff3 -m lite -t $prefix
    perl $LRSDAY_HOME/scripts/sort_gff3.pl -i $prefix.mitochondrial_genome.pre_sort.gff3 -t $prefix -o $prefix.mitochondrial_genome.gff3 -r $prefix.assembly.mitochondrial_genome.tidy.fa

    perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl -r $prefix.assembly.mitochondrial_genome.tidy.fa -g $prefix.mitochondrial_genome.gff3 -o $prefix.mitochondrial_genome.cds.fa
    perl $LRSDAY_HOME/scripts/cds2protein.pl -i $prefix.mitochondrial_genome.cds.fa -t $genetic_code_table -p $prefix.mitochondrial_genome

    perl $LRSDAY_HOME/scripts/prepare_PoFFgff_simple.pl -i $prefix.mitochondrial_genome.gff3 -o $prefix.mitochondrial_genome.PoFF.gff 
    perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i $prefix.mitochondrial_genome.trimmed_cds.fa -o $prefix.mitochondrial_genome.PoFF.ffn 
    perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i $prefix.mitochondrial_genome.pep.fa -o $prefix.mitochondrial_genome.PoFF.faa

    # clean up intermediate files
    if [[ $debug == "no" ]]
    then
	rm -r $$.tmp
	rm -r $prefix.mitochondrial_genome.pre_sort.gff3
	rm -r $prefix.assembly.mitochondrial_genome.fa
    fi
fi

        """
