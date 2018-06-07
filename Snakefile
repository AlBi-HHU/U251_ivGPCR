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
    params: mem_kb=config["biomart"]
    script:
        "scripts/getGO.R"

rule getKEGG:
    output:
        "KEGG/KEGG_pathways.csv"
    conda:
        "envs/python.yaml"
    params: mem_kb=config["getKEGG"]
    script:
        "scripts/getKEGG.py"


rule all_go_uniq:
# Make a list of unique GO terms
    input:
        "GO/GO_biomart.txt"
    output:
        "GO/all_go_uniq.txt"
    params: mem_kb=config["all_go_uniq"]
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
    params: mem_kb=config["ancestors"]
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
    params: mem_kb=config["topGO"]
    shell:
        "python scripts/mappingTopGO.py {input} > {output}"


rule counts:
    input:
        "CountTable_RNA_seq_U251_ivGPCR.xlsx"
    output:
        expand("count_files/{experiment}.txt", experiment=config["experiments"])
    conda:
      "envs/python.yaml"
    params: mem_kb=config["counts"]
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
    params: mem_kb=config["deseq2"]
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
    params: mem_kb=config["fit_BUM"]
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
    #conda:
        #"envs/python.yaml"
    params: mem_kb=config["run_heinz"]
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
    params: mem_kb=config["merge"]
    shell:
        "python scripts/merge.py {input[0]} {input[1]} GO/GO_biomart.txt GO/all_go_uniq_ancestors.txt > {output}"

rule module_to_study:
    input:
        "scores/{experiment}_{network}_{FDR}_module.mod"
    output:
        temp("{experiment}_{network}_{FDR}_modulestudy.txt")
    params: mem_kb=config["networktotxt"]
    shell:
        "python scripts/modtonames.py {input} {output}"

rule HPRDtotxt:
    input:
        "networks/HPRD_Release9_062910_Heinz.txt",
    output:
        "networks/HPRD_network.txt",
    params: mem_kb=config["networktotxt"]
    shell:
        r"cat {input} | sed  -e 's/\t/\n/g' | sort -u > {output}"

rule IREFKEGGtotxt:
    input:
        "networks/irefindex14+kegg09.txt"
    output:
        "networks/IREFKEGG_network.txt"
    params: mem_kb=config["networktotxt"]
    shell:
        r"cat {input} | sed  -e 's/\t/\n/g' | sort -u > {output}"

rule STRINGtotxt:
    input:
        "networks/9606.protein.links.v10_mapped_geneIDs_700.txt"
    output:
        "networks/STRING_network.txt"
    params: mem_kb=config["networktotxt"]
    shell:
        r"cat {input} | sed -e 's/\t/\n/g' | sort -u > {output}"

rule format_map:
    input:
        "GO/{experiment}_gsymb2go.map"
    output:
        "GO/{experiment}_gsymb2go_goaformat.map"
    params: mem_kb=config["format_map"]
    shell:
        "sed 's/,/;/g' {input} > {output}"

rule goa_to_eXamine:
    input:
        "scores/{experiment}_{network}_{FDR}_enrich.txt"
    output:
        "scores/{experiment}_{network}_{FDR}_sets.exm"
    params: mem_kb=config["goa_to_eXamine"]
    shell:
        "python scripts/enrich_to_eXamine.py {input} > {output}"
rule goaenrich:
    input:
        "GO/{experiment}_gsymb2go_goaformat.map",
        "{experiment}_{network}_{FDR}_modulestudy.txt",
        "networks/STRING_network.txt",
        "networks/HPRD_network.txt",
        "networks/IREFKEGG_network.txt"
    output:
        "scores/{experiment}_{network}_{FDR}_enrich.txt"
    conda:
        "envs/python.yaml"
    params: mem_kb=config["enrich"]
    shell:
        "python scripts/find_enrichment.py {input[1]} networks/{wildcards.network}_network.txt {input[0]} > {output}"

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
    params: mem_kb=config["eXamine_nodes"]
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
    params: mem_kb=config["eXamine_interactions"]
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
    params: mem_kb=config["eXamine_sets"]
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
    params: mem_kb=config["merge_heinz_deseq2"]
    shell:
        "touch {output}; scripts/merge_mod_deseq.py -m {input[0]} -d {input[1]} -o {output}"

rule filter_STRING:
    input:
        "networks/9606.protein.links.v10_mapped_geneIDs_700.txt"
    output:
        "networks/9606.protein.links.v10_mapped_geneIDs_700_no_UBC.txt"
    params: mem_kb=config["filters"]
    shell:
        "sed '/UBC/d' {input} > {output}"

rule filter_IREF:
    input:
        "networks/irefindex14+kegg09.txt"
    output:
        "networks/irefindex14+kegg09.txt_no_UBC.txt"
    params: mem_kb=config["filters"]
    shell:
        "sed '/UBC/d' {input} > {output}"

rule plot_stripchart:
    input:
        normcounts="deseq2/{experiment}_deseq2_normalized_counts.txt"
    output:
        "plots/{experiment}_top{n}_stripchart.pdf"
    conda:
        "envs/python.yaml"
    params: mem_kb=config["plot_stripchart"]
    script:
        "scripts/plot-stripchart.py"
