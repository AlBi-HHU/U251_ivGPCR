import pandas as pd
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("paper")
sns.set_style("ticks", {"grid.linestyle": ":"})


def plot_stripchart(counts, condition_re):
    counts = counts.stack().to_frame()
    counts = counts.reset_index()
    counts.columns = ["gene", "sample", "count"]
    condition = counts["sample"].str.extract(condition_re, expand=True)
    try:
        counts["condition"] = condition["condition"]
    except KeyError:
        raise ValueError("regular expression for condition does not contain a group 'condition'")
    if counts["condition"].isnull().any():
        raise ValueError("regular expression for condition did not match for sample {}".format(
            counts[~counts["condition"].isnull(), "gene"]))


    sns.stripplot(x="count", y="gene", hue="condition", data=counts, jitter=True, alpha=0.8)
    plt.xlabel("normalized counts")
    plt.grid(axis="y")
    plt.gca().set_xscale("log")



if __name__ == "__main__":
    counts = pd.read_table(snakemake.input.normcounts)
    n = int(snakemake.wildcards.n)
    plot_stripchart(counts.iloc[:n], snakemake.config["condition_re"])
    plt.savefig(snakemake.output[0], bbox_inches="tight")
