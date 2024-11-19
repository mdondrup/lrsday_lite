rule TE_annotation:
    conda: "../envs/gene_annotation.yaml"
    input: assembly="07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa",
           reannotate=".reannotate_setup"
    output:
        gff3="11.TE_Annotation/"+config["prefix"]+".nuclear_genome.TE.gff3",

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
debug={config[debug]}
threads={threads}        
        
source tools/reannotate.env        

clustalw=$CONDA_PREFIX/bin/clustalw2        
        
#######################################
# process the pipeline

## define a safe grep that will not result in pipefail if pattern not found        
c1grep() {{ grep "$@" || test $? = 1; }}

        
cd 11.TE_Annotation/        
        
TE_lib="$LRSDAY_HOME/data/TE_lib.v20221221.tidy.fa"
TE_lib_LTRonly="$LRSDAY_HOME/data/TE_lib.v20221221.LTRonly.tidy.fa"

echo $chrMT_tag | sed -e "s/;/\n/g" > $prefix.assembly.chrMT.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $genome_assembly -l $prefix.assembly.chrMT.list -m reverse -o $prefix.assembly.nuclear_genome.fa
perl $LRSDAY_HOME/scripts/tidy_fasta.pl -i $prefix.assembly.nuclear_genome.fa -o $prefix.assembly.nuclear_genome.tidy.fa
ln -fs $prefix.assembly.nuclear_genome.tidy.fa  $prefix.genome.fa

RepeatMasker -s -pa $threads -lib $TE_lib  $prefix.genome.fa  -xsmall  -gff
perl $reannotate_dir/REannotate_longname  -g -f $LRSDAY_HOME/data/fuzzy_defragmentation.txt  -k $clustalw -r 0.48764  -d 10000  -t $prefix.genome.fa.out  $prefix.genome.fa
rm -rf ${{prefix}}_REannotate_out        
mv -n REannotate_output ${{prefix}}_REannotate_out

cd ${{prefix}}_REannotate_out 

perl $LRSDAY_HOME/scripts/parse_REannotate_gff.pl -o $prefix.REannotate.raw.gff
mv $prefix.REannotate.raw.gff $prefix.REannotate.gff
perl $LRSDAY_HOME/scripts/parse_REannotate_out.pl -i $prefix.genome.fa.REannotation -p $prefix.TY_REannotate -r ./../$prefix.genome.fa

cat $prefix.TY_REannotate.complete.raw.fa  |c1grep ">" |c1grep  '(TY1|TY2)' |sed "s/>//g" >$prefix.TY_REannotate.complete.TY1TY2.raw.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.TY_REannotate.complete.raw.fa -l $prefix.TY_REannotate.complete.TY1TY2.raw.list  -o $prefix.TY_REannotate.complete.TY1TY2.raw.fa

cat $prefix.TY_REannotate.truncated.raw.fa  |c1grep ">" |c1grep  '(TY1|TY2)' |sed "s/>//g" >$prefix.TY_REannotate.truncated.TY1TY2.raw.list
perl $LRSDAY_HOME/scripts/select_fasta_by_list.pl -i $prefix.TY_REannotate.truncated.raw.fa -l $prefix.TY_REannotate.truncated.TY1TY2.raw.list  -o $prefix.TY_REannotate.truncated.TY1TY2.raw.fa

        
## The original code fails if grep doesn't find anything
## See: https://stackoverflow.com/a/49627999/326156        

        
for i in TY3 TY4 TY5 TSU4
  do
    echo $i -1   
    cat $prefix.TY_REannotate.complete.raw.fa  |  c1grep ">"  | c1grep "$i" | sed -e "s/>//g"  > $prefix.TY_REannotate.complete.$i.final.list
    echo $i - 2    
    cat $prefix.TY_REannotate.truncated.raw.fa  | c1grep ">" | c1grep "$i" | sed -e "s/>//g"  > $prefix.TY_REannotate.truncated.$i.final.list
  done

        
TY2_query="$LRSDAY_HOME/data/TY2_specific_region.fa"


db="$prefix.TY_REannotate.complete.TY1TY2.raw.fa"
db_tag="complete_TY1TY2_db"
if [ ! -s "$db" ]
then
    echo "$db is empty, skip .."
    touch $prefix.TY_REannotate.complete.TY1.final.list
    touch $prefix.TY_REannotate.complete.TY2.final.list
else
    makeblastdb  -in $db  -dbtype nucl -title $db_tag -hash_index -out $db_tag
    blastn -task blastn -query $TY2_query -num_threads $threads -db $db_tag -outfmt 7 >$prefix.$db_tag.blastn.fmt7.out
    perl $LRSDAY_HOME/scripts/filter_blast_result.pl -i $prefix.$db_tag.blastn.fmt7.out -pct_identity_cutoff 95 -aln_length_cutoff 500 -o $prefix.$db_tag.blastn.fmt7.I95L500.out
    cat $prefix.$db_tag.blastn.fmt7.I95L500.out |c1grep -v "^#" |cut -f 2  > $prefix.TY_REannotate.complete.TY2.list
    comm -23 <(sort $prefix.TY_REannotate.complete.TY1TY2.raw.list) <(sort $prefix.TY_REannotate.complete.TY2.list) > $prefix.TY_REannotate.complete.TY1.list
    cat $prefix.TY_REannotate.complete.TY1.list | sed "s/TY2/TY1/g" >$prefix.TY_REannotate.complete.TY1.final.list
    cat $prefix.TY_REannotate.complete.TY2.list | sed "s/TY1/TY2/g" >$prefix.TY_REannotate.complete.TY2.final.list
fi

db="$prefix.TY_REannotate.truncated.TY1TY2.raw.fa"
db_tag="truncated_TY1TY2_db"
if [ ! -s "$db" ]
then
    echo "$db is empty, skip .."  
    touch $prefix.TY_REannotate.truncated.TY1.final.list
    touch $prefix.TY_REannotate.truncated.TY2.final.list
else
    makeblastdb  -in $db  -dbtype nucl -title $db_tag -hash_index -out $db_tag
    blastn -task blastn -query $TY2_query -num_threads $threads -db $db_tag -outfmt 7 >$prefix.$db_tag.blastn.fmt7.out
    perl $LRSDAY_HOME/scripts/filter_blast_result.pl -i $prefix.$db_tag.blastn.fmt7.out -pct_identity_cutoff 95 -aln_length_cutoff 100 -o $prefix.$db_tag.blastn.fmt7.I95L100.out
    cat $prefix.$db_tag.blastn.fmt7.I95L100.out | c1grep -v "^#" |cut -f 2 > $prefix.TY_REannotate.truncated.TY2.list
    comm -23 <(sort $prefix.TY_REannotate.truncated.TY1TY2.raw.list) <(sort $prefix.TY_REannotate.truncated.TY2.list) > $prefix.TY_REannotate.truncated.TY1.list
    cat $prefix.TY_REannotate.truncated.TY1.list |sed "s/TY2/TY1/g" >$prefix.TY_REannotate.truncated.TY1.final.list
    cat $prefix.TY_REannotate.truncated.TY2.list |sed "s/TY1/TY2/g" >$prefix.TY_REannotate.truncated.TY2.final.list
fi

LTR_query="$prefix.TY_REannotate.soloLTR.raw.fa"
db=$TE_lib_LTRonly
db_tag="soloLTR_db";
if [ ! -s "$LTR_query" ]
then
    echo "$LTR_query is empty, skip .."   
    touch $prefix.TY_soloLTR.refined.nr.gff
else
    makeblastdb -in $db -dbtype nucl -title $db_tag -hash_index -out $db_tag
    blastn -task blastn -query $LTR_query -num_threads $threads -db $db_tag -outfmt 7 >$prefix.$db_tag.soloLTR.blastn.fmt7.out
    perl $LRSDAY_HOME/scripts/trim_soloLTR_by_blast.pl -q $LTR_query -b  $prefix.$db_tag.soloLTR.blastn.fmt7.out -p $prefix  -i 75 -l 100
    bedtools sort  -i  $prefix.TY_REannotate.soloLTR.refined.gff > $prefix.TY_REannotate.soloLTR.refined.sorted.gff
    perl $LRSDAY_HOME/scripts/rm_overlap_features_from_gff_simple.pl  -r ./../$prefix.genome.fa  -i $prefix.TY_REannotate.soloLTR.refined.sorted.gff -o $prefix.TY_soloLTR.refined.nr.gff 
fi

for i in TY1 TY2 TY3 TY4 TY5 TSU4
do
    perl $LRSDAY_HOME/scripts/TY_list2gff3.pl -i $prefix.TY_REannotate.complete.$i.final.list -o $prefix.TY_REannotate.complete.$i.final.gff -t $prefix
        
    perl $LRSDAY_HOME/scripts/TY_list2gff3.pl -i $prefix.TY_REannotate.truncated.$i.final.list -o $prefix.TY_REannotate.truncated.$i.final.gff -t $prefix
done

cat $prefix.TY_REannotate.*.final.gff > $prefix.TY.complete_plus_truncated.final.gff

bedtools sort -i $prefix.TY.complete_plus_truncated.final.gff > $prefix.TY.complete_plus_truncated.final.sorted.gff 
bedtools intersect -v -a $prefix.TY_soloLTR.refined.nr.gff -b $prefix.TY.complete_plus_truncated.final.sorted.gff  >$prefix.TY.soloLTR.final.gff

cat $prefix.TY.complete_plus_truncated.final.gff $prefix.TY.soloLTR.final.gff > $prefix.TY.all.final.gff

perl $LRSDAY_HOME/scripts/tidy_TE_gff3.pl -r ./../$prefix.genome.fa -i $prefix.TY.all.final.gff -o $prefix.nuclear_genome.TE.raw.gff3 -t $prefix 
perl $LRSDAY_HOME/scripts/adjust_TY_annotation_for_gff3.pl -i $prefix.nuclear_genome.TE.raw.gff3 -o ./../$prefix.nuclear_genome.TE.gff3 


if [[ $debug == "no" ]]
then
    rm -rf $prefix.*.final.gff
    rm -rf $prefix.nuclear_genome.TE.raw.gff3
fi

cd ..

# clean up intermediate files
if [[ $debug == "no" ]]
then
    rm -rf $prefix.assembly.nuclear_genome*
    rm -rf $prefix.genome.fa*
fi
        
        """
                 

             
        
        
    
