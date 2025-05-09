#### Alternative TE annotation with RepeatMasker using the sequences from:
####
#### Czaja W, Bensasson D, Ahn HW, Garfinkel DJ, Bergman CM (2020)
#### Evolution of Ty1 copy number control in yeast by horizontal transfer and recombination.
#### PLOS Genetics 16(2): e1008632. https://doi.org/10.1371/journal.pgen.1008632

#### WU-Blast can be replaced with RPM-Blast or AB-Blast

db_url=config["alternative_te_url"] # configure the URL
engine=config["alt_rm_blast_engine"] or "rmblast"
### run only if the alternative url is configured
if (db_url):    
    rule download_db:
        output: "11A.Alternative_TE_Annotation/TY_DB.fna"
        params: dburl=db_url
        log: "logs/++download_TE.log"        
        shell:
            r"""
            mkdir -p 11A.Alternative_TE_Annotation/
            cd 11A.Alternative_TE_Annotation/
            wget -c {params.dburl} -O ./TY_DB.fna -o ../{log}
            """

    rule alt_repmasker:
        conda: "../envs/gene_annotation.yaml"
        input:
            #rm = ".repmasker_setup",
            #wu = ".wublast_rm_setup"  if (engine == "wublast") else [],
            lib = "11A.Alternative_TE_Annotation/TY_DB.fna",
            asm = "07.Supervised_Final_Assembly/"+ config["prefix"] +".assembly.final.fa"
        output:
            xm="11A.Alternative_TE_Annotation/"+config["prefix"]+".assembly.final.fa.out.xm",
            out="11A.Alternative_TE_Annotation/"+config["prefix"]+".assembly.final.fa.out",
            masked="11A.Alternative_TE_Annotation/"+config["prefix"]+".assembly.final.fa.masked",
            gff="11A.Alternative_TE_Annotation/"+config["prefix"]+".assembly.final.fa.out.gff"
        threads: 120
        params:
            engine=engine
        log:
            "logs/" + config["prefix"] + ".repeatmasker.log"
        message:
            """
            Running alternative TE annotation with RepeatMasker using {params.engine} engine. 
            """
        shell:
            """
            RepeatMasker -lib {input.lib} -s -xsmall -nolow -no_is -pa {threads} -gff -html -xm -e {params.engine} -dir $(dirname {output.out}) {input.asm} > {log} 2>&1
            """

    rule tyson:
        conda: "../envs/gene_annotation.yaml"
        input:
            rmout="11A.Alternative_TE_Annotation/"+config["prefix"]+".assembly.final.fa.out",
            masked="11A.Alternative_TE_Annotation/"+config["prefix"]+".assembly.final.fa.masked"
        output:
            gff="11A.Alternative_TE_Annotation/"+config["prefix"]+".tyson.gff",
            bed="11A.Alternative_TE_Annotation/"+config["prefix"]+".tyson.bed",
            ty1bed="11A.Alternative_TE_Annotation/"+config["prefix"]+".TY1complete.bed",
            counts="11A.Alternative_TE_Annotation/"+config["prefix"]+".counts",
            fa="11A.Alternative_TE_Annotation/"+config["prefix"]+".TyAll.fa",
            ty1fa="11A.Alternative_TE_Annotation/"+config["prefix"]+".TY1complete.fa",
        params:
            prefix=config["prefix"]
        log:
            "logs/" + config["prefix"] + ".tyson.log"
        message:
            r"""
            Running Ty element classification similar to Czaja et al. (2020).
            
            """ 
        shell:
            r"""
            cd 11A.Alternative_TE_Annotation/
            P='{params.prefix}_{{nr}}'
            perl ../workflow/scripts/TySon.pl $(basename {input.rmout}) --no-bed --counts --prefix {params.prefix} > $(basename {output.gff}) 2> ../{log}
             perl ../workflow/scripts/TySon.pl $(basename {input.rmout}) --bed --no-counts --prefix {params.prefix} > $(basename {output.bed}) 2>> ../{log}
            bedtools getfasta -s -fi $(basename {input.masked}) -bed $(basename {output.bed}) -nameOnly  \
             >  $(basename {output.fa}) 2>> ../{log}
            perl ../workflow/scripts/TySon.pl $(basename {input.rmout}) --bed --class_filter 'TY1' --complete_only --no-counts --prefix {params.prefix} > $(basename {output.ty1bed}) 2>> ../{log}
            bedtools getfasta -s -fi $(basename {input.masked}) -bed $(basename {output.ty1bed}) -nameOnly  \
             >  $(basename {output.ty1fa}) 2>> ../{log}
            

            """
            
