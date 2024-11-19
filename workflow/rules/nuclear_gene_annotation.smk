### Nuclear genome annotation

rule nuclear_gene_annotation:    
    input: assembly="07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa",
           repmsker=".repmasker_setup"
    output:
        EVMDIR= directory("09.Nuclear_Gene_Annotation/"+config["prefix"]+".nuclear_genome.EVM.output"), # make sure this is removed on error 
        gff3= "09.Nuclear_Gene_Annotation/"+config["prefix"]+".nuclear_genome.gff3",
        cdsfa= "09.Nuclear_Gene_Annotation/"+config["prefix"]+".nuclear_genome.cds.fa",
         pepfa= "09.Nuclear_Gene_Annotation/"+config["prefix"]+".nuclear_genome.pep.fa",
         pofffaa="09.Nuclear_Gene_Annotation/"+config["prefix"]+".nuclear_genome.PoFF.faa",
         poffgff="09.Nuclear_Gene_Annotation/"+config["prefix"]+".nuclear_genome.PoFF.gff"
    conda: "../envs/gene_annotation.yaml"
    threads: 60
    shell:
        r"""

#######################################
# set project-specific variables
#        
export LRSDAY_HOME=$(realpath LRSDAY)
prefix={config[prefix]} # The file name prefix (only allowing strings of alphabetical letters, numbers, and underscores) for the processing sample. Default = "CPG_1a" for the testing example.
genome_assembly=$(realpath {input.assembly}) # The path of the input genome assembly.
chrMT_tag={config[chrMT]} # The sequence name for the mitochondrial genome in the final assembly. If there are multiple sequences, use a single ';' to separate them. e.g. "chrMT_part1;chrMT_part2". Default = "chrMT".

export EVM_HOME=$(realpath $CONDA_PREFIX/opt/evidencemodeler-1.*/)
export bin_dir=$(dirname $(which maker))        
export LIBDIR=$CONDA_PREFIX/share/RepeatMasker/Libraries # RepeatMasker Database
        
threads={threads} # The number of threads to use. Default = "8".
maker_opts="$LRSDAY_HOME/misc/maker_opts.customized.ctl" # The configuration file for MAKER. You can edit this file if you have native transciptome/EST data for the strain/species that you sequenced or if you want to adapt it to annotate other eukaryotic organisms. Otherwise, please keep it unchanged. Please note that if this file is in the same directory where this bash script is executed, the file name cannot be "maker_opts.ctl".
EVM_weights="$LRSDAY_HOME/misc/EVM_weights.customized.txt" # The configuration file for EVM. A list of numeric weight values to be applied to each type of evidence.
debug={config[debug]} # Whether to keep intermediate files for debugging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".bioperl_exonerate.yaml

#######################################
# process the pipeline

cd 09.Nuclear_Gene_Annotation/
        
genome_tag="$prefix"

echo "genome_tag=$genome_tag"
echo "genome_assembly=$genome_assembly"

# convert the genome assembly file to all uppercases

echo $chrMT_tag | sed -e "s/;/\n/g" > $genome_tag.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $genome_tag.assembly.chrMT.list -m reverse -o $genome_tag.assembly.nuclear_genome.fa
# perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $genome_tag.assembly.chrMT.list -m normal -o $genome_tag.assembly.mitochondrial_genome.fa
perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $genome_tag.assembly.nuclear_genome.fa -o $genome_tag.assembly.nuclear_genome.tidy.fa

cp $maker_opts maker_opts.ctl
#cp $LRSDAY_HOME/misc/maker_exe.ctl . ## All tools in PATH
cp ../config/maker_exe.ctl .
cp $LRSDAY_HOME/misc/maker_bopts.ctl .
cp $LRSDAY_HOME/misc/maker_evm.ctl .

maker -fix_nucleotides -genome $genome_tag.assembly.nuclear_genome.tidy.fa -cpus $threads -base $genome_tag

# $maker_dir/fasta_merge -d $genome_tag.maker.output/${{genome_tag}}_master_datastore_index.log -o $genome_tag.nuclear_genome.maker.fasta
gff3_merge -d $genome_tag.maker.output/${{genome_tag}}_master_datastore_index.log -n -g -o $genome_tag.nuclear_genome.maker.raw.gff3
perl $LRSDAY_HOME/scripts/collect_maker_evidences.pl -t $genome_tag -p $genome_tag.nuclear_genome
cat $genome_tag.nuclear_genome.maker.raw.gff3 |egrep "=trnascan" > $genome_tag.nuclear_genome.maker.raw.tRNA.gff3
# cat $genome_tag.nuclear_genome.maker.raw.gff3 |egrep "=snoscan" > $genome_tag.nuclear_genome.maker.raw.snoscan.gff3 # developmental feature
cat $genome_tag.nuclear_genome.maker.raw.gff3 |egrep -v "=trnascan"|egrep -v "=snoscan" > $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.gff3

# use EVM to further polishing the annotation for multi-exon genes
perl $LRSDAY_HOME/scripts/filter_gff3_for_single_exon_genes.pl -i $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.gff3 -o $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.single_exon_gene.gff3

bedtools intersect -v  -a $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.gff3 \
		       -b $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.single_exon_gene.gff3 \
		       > $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.multiple_exon_gene.gff3
bedtools intersect -v  -a $genome_tag.nuclear_genome.protein_evidence.gff3 \
		       -b $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.single_exon_gene.gff3 \
		       > $genome_tag.nuclear_genome.protein_evidence.for_gene_model_refinement.gff3
bedtools intersect -v  -a $genome_tag.nuclear_genome.est_evidence.gff3 \
		       -b $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.single_exon_gene.gff3 \
		       > $genome_tag.nuclear_genome.est_evidence.for_gene_model_refinement.gff3

perl $LRSDAY_HOME/scripts/exonerate2gene_maker.pl -i $genome_tag.nuclear_genome.protein_evidence.for_gene_model_refinement.gff3 -o $genome_tag.nuclear_genome.protein_evidence.complementary_gene_model.gff3 
perl $LRSDAY_HOME/scripts/exonerate2gene_maker.pl -i $genome_tag.nuclear_genome.est_evidence.for_gene_model_refinement.gff3 -o $genome_tag.nuclear_genome.est_evidence.complementary_gene_model.gff3

cat $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.gff3 $genome_tag.nuclear_genome.protein_evidence.complementary_gene_model.gff3 $genome_tag.nuclear_genome.est_evidence.complementary_gene_model.gff3 > $genome_tag.nuclear_genome.maker.combined.gff3
set -x
mkdir -p $genome_tag.nuclear_genome.EVM.output
cd $genome_tag.nuclear_genome.EVM.output
$EVM_HOME/EvmUtils/partition_EVM_inputs.pl \
    --genome ./../$genome_tag.assembly.nuclear_genome.tidy.fa \
    --gene_predictions ./../$genome_tag.nuclear_genome.maker.combined.gff3 \
    --protein_alignments ./../$genome_tag.nuclear_genome.protein_evidence.gff3 \
    --transcript_alignments ./../$genome_tag.nuclear_genome.est_evidence.gff3 \
    --segmentSize 100000 \
    --overlapSize 10000 \
    --partition_listing $genome_tag.nuclear_genome.partitions_list.out

$EVM_HOME/EvmUtils/write_EVM_commands.pl \
    --genome ./../$genome_tag.assembly.nuclear_genome.tidy.fa \
    --weights $EVM_weights  \
    --gene_predictions ./../$genome_tag.nuclear_genome.maker.combined.gff3 \
    --output_file_name $genome_tag.nuclear_genome.evm.out \
    --partitions $genome_tag.nuclear_genome.partitions_list.out \
    > $genome_tag.nuclear_genome.commands.list

$EVM_HOME/EvmUtils/execute_EVM_commands.pl $genome_tag.nuclear_genome.commands.list | tee $genome_tag.nuclear_genome.EVM_run.log
$EVM_HOME/EvmUtils/recombine_EVM_partial_outputs.pl \
    --partitions $genome_tag.nuclear_genome.partitions_list.out \
    --output_file_name $genome_tag.nuclear_genome.evm.out
$EVM_HOME/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
    --partitions $genome_tag.nuclear_genome.partitions_list.out \
    --output $genome_tag.nuclear_genome.evm.out \
    --genome ./../$genome_tag.assembly.nuclear_genome.tidy.fa

perl $LRSDAY_HOME/scripts/collect_EVM_gff3.pl -p $genome_tag.nuclear_genome  -r ./../$genome_tag.assembly.nuclear_genome.tidy.fa

cat $genome_tag.nuclear_genome.EVM.raw.gff3 ./../$genome_tag.nuclear_genome.maker.raw.tRNA.gff3 > $genome_tag.nuclear_genome.EVM.raw.with_tRNA.gff3
perl $LRSDAY_HOME/scripts/tidy_maker_gff3.pl \
    -i $genome_tag.nuclear_genome.EVM.raw.with_tRNA.gff3 \
    -r ./../$genome_tag.assembly.nuclear_genome.tidy.fa \
    -t $genome_tag \
    -o $genome_tag.nuclear_genome.EVM.gff3
cp $genome_tag.nuclear_genome.EVM.gff3 ./../$genome_tag.nuclear_genome.gff3
cd ..

perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl \
    -r $genome_tag.assembly.nuclear_genome.tidy.fa \
    -g $genome_tag.nuclear_genome.gff3 \
    -o $genome_tag.nuclear_genome.cds.fa
perl $LRSDAY_HOME/scripts/cds2protein.pl \
    -i $genome_tag.nuclear_genome.cds.fa \
    -p $genome_tag.nuclear_genome \
    -t 1


# filtered out snoRNA annotation since it is still an experimental features suffering from redundant annotations
cp $genome_tag.nuclear_genome.gff3 $genome_tag.nuclear_genome.gff3.tmp
cat $genome_tag.nuclear_genome.gff3.tmp | egrep -v "snoRNA" > $genome_tag.nuclear_genome.snoRNA_filtered.gff3
perl $LRSDAY_HOME/scripts/label_pseudogene_in_gff3.pl -i $genome_tag.nuclear_genome.snoRNA_filtered.gff3 -l $genome_tag.nuclear_genome.manual_check.list -o $genome_tag.nuclear_genome.gff3

perl $LRSDAY_HOME/scripts/extract_cds_from_tidy_gff3.pl \
    -r $genome_tag.assembly.nuclear_genome.tidy.fa \
    -g $genome_tag.nuclear_genome.gff3 \
    -o $genome_tag.nuclear_genome.cds.fa
perl $LRSDAY_HOME/scripts/cds2protein.pl \
    -i $genome_tag.nuclear_genome.cds.fa \
    -p $genome_tag.nuclear_genome \
    -t 1

perl $LRSDAY_HOME/scripts/prepare_PoFFgff_simple.pl -i $genome_tag.nuclear_genome.gff3 -o $genome_tag.nuclear_genome.PoFF.gff 
perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i $genome_tag.nuclear_genome.trimmed_cds.fa -o $genome_tag.nuclear_genome.PoFF.ffn 
perl $LRSDAY_HOME/scripts/prepare_PoFFfaa_simple.pl -i $genome_tag.nuclear_genome.pep.fa -o $genome_tag.nuclear_genome.PoFF.faa

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm $genome_tag.assembly.nuclear_genome.fa
    rm $genome_tag.nuclear_genome.maker.raw.tRNA.gff3
    rm $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.gff3
    rm $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.single_exon_gene.gff3
    rm $genome_tag.nuclear_genome.maker.raw.protein_coding_gene.multiple_exon_gene.gff3
    rm $genome_tag.nuclear_genome.protein_evidence.for_gene_model_refinement.gff3
    rm $genome_tag.nuclear_genome.est_evidence.for_gene_model_refinement.gff3
    rm $genome_tag.nuclear_genome.protein_evidence.complementary_gene_model.gff3
    rm $genome_tag.nuclear_genome.est_evidence.complementary_gene_model.gff3
    rm $genome_tag.nuclear_genome.gff3.tmp
    rm $genome_tag.nuclear_genome.maker.combined.gff3
    rm $genome_tag.nuclear_genome.snoRNA_filtered.gff3
    rm -rf _Inline
    rm *.ctl
fi

        """        
