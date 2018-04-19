configfile: "config.yaml" # json oder yaml file

targets = ["scores/{}_{}_{}_module_fc_pval.txt".format(experiment, network, fdr) for experiment in config["experiments"] for network in config["networks"] for fdr in config["FDRs"][network][experiment]] + ["GO/{}_gsymb2go.map".format(experiment) for experiment in config["experiments"]] + ["data-sets/{}_{}_{}/interactions.links".format(experiment, network, fdr) for experiment in config["experiments"] for network in config["networks"] for fdr in config["FDRs"][network][experiment]] + ["data-sets/{}_{}_{}/go_and_kegg.memberships".format(experiment, network, fdr) for experiment in config["experiments"] for network in config["networks"] for fdr in config["FDRs"][network][experiment]]

target_modules = ["scores/{}_{}_{}_module_{}.txt".format(experiment, network, fdr, i) for experiment in config["experiments"] for network in config["networks"] for fdr in config["FDRs"][network][experiment] for i in range(1) ]

target_stripcharts = expand("plots/{experiment}_top{n}_stripchart.pdf", experiment=config["experiments"], n=20)


#debug
#print(target_modules)

wildcard_constraints:
    experiment = '\w+_vs_\w+'

rule all:
    input:
        targets,
#        target_modules,
        "GO/all_go_uniq_ancestors.txt", #, "GO/all_go_uniq_ancestors.txt"
        "KEGG/KEGG_pathways.csv",
        target_stripcharts


rule biomart:
    output:
        "GO/GO_biomart.txt"
    conda:
        "envs/go.yaml"
    #params:
    #    version=config["ENSEMBL"]["version"]
    script:
        "scripts/getGO.R"

rule getKEGG:
    output:
        "KEGG/KEGG_pathways.csv"
    conda:
        "envs/python.yaml"
    script:
        "scripts/getKEGG.py"


rule all_go_uniq:
# Make a list of unique GO terms
    input:
        "GO/GO_biomart.txt"
    output:
        "GO/all_go_uniq.txt"
    shell:
        "cut -f2,4 {input} | sort -u  > {output}"

rule ancestors:
# extend the list of unique GO terms with their ancestors
    input:
        "GO/all_go_uniq.txt"
    output:
        "GO/all_go_uniq_ancestors.txt"
    conda:
        "envs/go.yaml"
    shell:
        "Rscript scripts/getAncestors.R {input} | sort -u > {output}"

rule topGO:
# create, for each experiment, a list of occuring genes with all the GO terms they occur in
    input:
        "GO/GO_biomart.txt",
        "scores/{experiment}_pvals.txt",
        "GO/all_go_uniq_ancestors.txt"
    output:
        "GO/{experiment}_gsymb2go.map"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/mappingTopGO.py {input} > {output}"


rule counts:
    input:
        "CountTable_RNA_seq_U251_ivGPCR.xlsx"
    output:
        expand("count_files/{experiment}.txt", experiment=config["experiments"])
    conda:
      "envs/python.yaml"
    script: "./scripts/xls_to_count_file.py"

rule deseq2:
    input:
      matrix="design_matrices/design_matrix_{experiment}.txt",
      counts="count_files/{experiment}.txt"
    output:
      pvals="deseq2/{experiment}_deseq2.txt",
      hist="deseq2/{experiment}_p_val_distr.pdf",
      normcounts="deseq2/{experiment}_deseq2_normalized_counts.txt"
    conda:
      "envs/deseq2.yaml"
    script:
        "scripts/DESeq2_DEA.R"

rule fit_BUM:
    input:
        deseq_out="deseq2/{experiment}_deseq2.txt"
    output:
        "scores/{experiment}_pvals.txt",
        "scores/{experiment}_bum_fit.txt"
    conda:
        "envs/bionet.yaml"
    script:
        "scripts/fitBUM.R"


# rule get_workingcopy_network:
#     input: network=snakemake.wildcards.network
#     output: "scores/{experiment}_{network}_{FDR}_network.txt"
#     shell:
#         "cp {input} {output}"

rule run_heinz:
    input:
        pvals="scores/{experiment}_pvals.txt",
        bum_params="scores/{experiment}_bum_fit.txt",
        #network = "scores/{experiment}_{network}_{FDR}_network.txt"
    output:
        "scores/{experiment}_{network}_{FDR}_module.txt",
        "scores/{experiment}_{network}_{FDR}_module.mod",
        "scores/{experiment}_{network}_{FDR}_module.pdf"
        # "scores/{experiment}_{network}_{FDR}_module_0.txt",
        # "scores/{experiment}_{network}_{FDR}_module_0.mod",
        # "scores/{experiment}_{network}_{FDR}_module_0.pdf"
    conda:
        "envs/python.yaml"
    script:
        "scripts/run_heinz.py"

rule merge:
    input:
        "scores/{experiment}_pvals.txt",
        full_module="scores/{experiment}_{network}_{FDR}_module.txt"
    output:
        "scores/{experiment}_{network}_{FDR}_nodes.exm"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/merge.py {input[0]} {input[1]} GO/GO_biomart.txt GO/all_go_uniq_ancestors.txt > {output}"

rule enrich:
    input: "GO/{experiment}_gsymb2go.map",
        "scores/{experiment}_{network}_{FDR}_module.mod"
    output: "scores/{experiment}_{network}_{FDR}_sets.exm",
        temp("tmp/{experiment}_{network}_{FDR}_enr_GO_P.txt"),
        temp("tmp/{experiment}_{network}_{FDR}_enr_GO_F.txt"),
        temp("tmp/{experiment}_{network}_{FDR}_enr_GO_C.txt")
    conda:
        "envs/topgo.yaml"
    shell:
        """
        Rscript scripts/enrichment.R {input[0]} {input[1]} BP {output[1]}
        Rscript scripts/enrichment.R {input[0]} {input[1]} MF {output[2]}
        Rscript scripts/enrichment.R {input[0]} {input[1]} CC {output[3]}
        cat {output[1]} > {output[0]}
        tail -n +2 {output[2]} >> {output[0]}
        tail -n +2 {output[3]} >> {output[0]}
        """

rule eXamine_nodes:
    input:
        "scores/{experiment}_{network}_{FDR}_module.mod",
        "scores/{experiment}_{network}_{FDR}_module_fc_pval.txt"
    output:
        "data-sets/{experiment}_{network}_{FDR}/modules.memberships",
        "data-sets/{experiment}_{network}_{FDR}/proteins.nodes",
        "data-sets/{experiment}_{network}_{FDR}/modules.annotations"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/proteins.py -m {input[0]} -fc {input[1]} -ol {output[0]} -on {output[1]} -om {output[2]}"

rule eXamine_interactions:
    input:
        lambda wildcards: os.path.join("networks", config["networks"][wildcards.network]),
        "data-sets/{experiment}_{network}_{FDR}/proteins.nodes"
    output:
        "data-sets/{experiment}_{network}_{FDR}/interactions.links"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/interactions.py -e {input[0]} -n {input[1]} -o {output}"

rule eXamine_sets:
    input:
        "GO/GO_biomart.txt",
        "data-sets/{experiment}_{network}_{FDR}/proteins.nodes",
        "scores/{experiment}_{network}_{FDR}_sets.exm",
        "KEGG/KEGG_pathways.csv",
        "scores/{experiment}_{network}_{FDR}_module.mod",
        "scores/{experiment}_{network}_{FDR}_module.txt"
    output:
        "data-sets/{experiment}_{network}_{FDR}/go_and_kegg.annotations",
        "data-sets/{experiment}_{network}_{FDR}/go_and_kegg.memberships"
    conda:
        "envs/python.yaml"
    shell:
        """
        python scripts/go_modules.py  -g {input[0]} -n {input[1]} -s {input[2]} -oa {output[0]} -ol {output[1]}
        python scripts/analyzeKEGG.py -m {input[4]} -b {input[5]} -ao {output[0]} -mo {output[1]}
        """

rule merge_heinz_deseq2:
    input:
        module = "scores/{experiment}_{network}_{FDR}_module.mod",
        deseq2_out = "deseq2/{experiment}_deseq2.txt"
    output:
        merge_res = "scores/{experiment}_{network}_{FDR}_module_fc_pval.txt"
    conda:
        "envs/python.yaml"
    shell:
        "touch {output}; scripts/merge_mod_deseq.py -m {input[0]} -d {input[1]} -o {output}"

rule filter_STRING:
    input:
        "networks/9606.protein.links.v10_mapped_geneIDs_700.txt"
    output:
        "networks/9606.protein.links.v10_mapped_geneIDs_700_no_UBC.txt"
    shell:
        "sed '/UBC/d' {input} > {output}"

rule filter_IREF:
    input:
        "networks/irefindex14+kegg09.txt"
    output:
        "networks/irefindex14+kegg09.txt_no_UBC.txt"
    shell:
        "sed '/UBC/d' {input} > {output}"

rule plot_stripchart:
    input:
        normcounts="deseq2/{experiment}_deseq2_normalized_counts.txt"
    output:
        "plots/{experiment}_top{n}_stripchart.pdf"
    conda:
        "envs/python.yaml"
    script:
        "scripts/plot-stripchart.py"
