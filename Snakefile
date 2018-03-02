configfile: "config.yaml" # json oder yaml file

targets = ["scores/{}_{}_{}_module_fc_pval.txt".format(experiment, fdr, network) for experiment in config["experiments"] for network in config["networks"] for fdr in config["FDRs"][network][experiment]] + ["GO/{}_gsymb2go.map".format(experiment) for experiment in config["experiments"]] + ["data-sets/{}_{}_{}/interactions.links".format(experiment, fdr, network) for experiment in config["experiments"] for network in config["networks"] for fdr in config["FDRs"][network][experiment]]


#debug
print(targets)

wildcard_constraints:
    experiment = '\w+_vs_\w+'

rule all:
    input: targets, "GO/all_go_uniq_ancestors.txt"#, "GO/all_go_uniq_ancestors.txt"

rule biomart:
    output:
        "GO/GO_biomart.txt"
    run:
        import time
        excep = True
        while excep:
            try:
                shell("Rscript scripts/getGO.R > {output}")
                excep = False
            except Exception:
                time.sleep(5)

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
    shell:
        "python scripts/mappingTopGO.py {input} > {output}"


rule counts:
    input: "CountTable_RNA_seq_U251_ivGPCR.xlsx"
    output: expand("count_files/{experiment}.txt", experiment=config["experiments"])
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

rule run_heinz:
    input:
        pvals="scores/{experiment}_pvals.txt",
        bum_params="scores/{experiment}_bum_fit.txt"
    output:
        "scores/{experiment}_{FDR}_{network}_module.txt",
        "scores/{experiment}_{FDR}_{network}_module.mod",
        "scores/{experiment}_{FDR}_{network}_module.pdf"
    script:
        "scripts/run_heinz.py"

rule heinz_res:
    input: "scores/{experiment}_{FDR}_{network}_module.txt"
    output: "scores/{experiment}_{FDR}_{network}_module.res"
    shell: "grep -v NaN {input} | grep -v \"#\" > {output}"

rule merge:
    input: "scores/{experiment}_pvals.txt",
        full_module="scores/{experiment}_{FDR}_{network}_module.txt"
    output: "scores/{experiment}_{FDR}_{network}_nodes.exm"
    shell:
        "python scripts/merge.py {input[0]} {input[1]} GO/GO_biomart.txt GO/all_go_uniq_ancestors.txt > {output}"

rule enrich:
    input: "GO/{experiment}_gsymb2go.map",
        "scores/{experiment}_{FDR}_{network}_module.res"
    output: "scores/{experiment}_{FDR}_{network}_sets.exm",
        temp("tmp/{experiment}_{FDR}_{network}_enr_GO_P.txt"),
        temp("tmp/{experiment}_{FDR}_{network}_enr_GO_F.txt"),
        temp("tmp/{experiment}_{FDR}_{network}_enr_GO_C.txt"),
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
    input: "scores/{experiment}_{FDR}_{network}_module.res"
    output: "data-sets/{experiment}_{FDR}_{network}/modules.links",
        "data-sets/{experiment}_{FDR}_{network}/proteins.nodes",
        "data-sets/{experiment}_{FDR}_{network}/modules.annotations"
    shell: "python scripts/proteins.py -m {input} -ol {output[0]} -on {output[1]} -om {output[2]}"

rule eXamine_interactions:
    input: lambda wildcards: os.path.join("networks", config["networks"][wildcards.network]),
    output: "data-sets/{experiment}_{FDR}_{network}/interactions.links"
    shell:
        "python scripts/interactions.py -e {network} -o {output}"


rule merge_heinz_deseq2:
    input:
        module = "scores/{experiment}_{FDR}_{network}_module.mod",
        deseq2_out = "deseq2/{experiment}_deseq2.txt"
    output:
        merge_res = "scores/{experiment}_{FDR}_{network}_module_fc_pval.txt"
    shell:
        "touch {output}; scripts/merge_mod_deseq.py -m {input[0]} -d {input[1]} -o {output}"

rule filter_STRING:
    input:
        "networks/9606.protein.links.v10_mapped_geneIDs_700.txt"
    output:
        "networks/9606.protein.links.v10_mapped_geneIDs_700_no_UBC.txt"
    shell:
        "sed '/UBC/d' {input} > {output}"

# rule getGO:
#     input:
#     output: "GO/GO-biomart.txt"
#     conda: "envs/GO.yaml"
#     script:
#         "scripts/getGO.R"
#
# rule getGOancestors:
#     input: "GO/GO-biomart.txt"
#     output: "GO/all_go_uniq_ancestors.txt"
#     conda: "envs/GO.yaml"
#     shell: "Rscript scripts/getAncestors.R {input} > {output}"
#
# rule topGO:
#     input:
#         "GO/GO-biomart.txt",
#         "scores/{experiment}_pvals.txt",
#         "GO/all_go_uniq_ancestors.txt"
#     output:
#         "GO/{experiment}_gsymb2go.map"
#     script:
#         "scripts/mappingTopGO.py"
#
# rule merge:
#     input:
#         "scores/{experiment}_pvals.txt",
#         "scores/{experiment}_{FDR}_{network}_module.txt",
#         #"GO/GO-biomart.txt",
#         #"GO/all_go_uniq_ancestors.txt"
#     output:
#         "eXamine/{experiment}_{FDR}_{network}_nodes.txt"
#     script:
#         "scripts/merge.py"

#         # <INPUT_FILE_PVAL> <INPUT_FILE_MODULE> <INPUT_GO_FILE> <INPUT_GO_ANCESTOR> <OUTPUT_FILE_NAME>
#
# # rule merge:
# #     input:
# #         "scores/stimulated_vs_US28_pvals.txt",
# #         "scores/stimulated_vs_US28_1e-175_STRING_module.txt",
# #         "GO/GO-biomart.txt",
# #         "GO/all_go_uniq_ancestors.txt"
# #     output:
# #         "eXamine/stimulated_vs_US28_1e-175_STRING_nodes.txt"
# #     script:
# #         "scripts/merge.py"
