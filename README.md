# ESPIDAM Tutorial 2026

This repository contains a three-part tutorial on **network-based epidemic modelling**, built around an agent-based SIR/SIHR framework in Julia. You will progressively layer realism onto the model:

| Notebook         | Theme                                                                                                                                                                                          |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Part 1** | Synthetic networks (random, small-world, preferential attachment) + basic SIR                                                                                                                  |
| **Part 2** | Empirical degree distributions (configuration model, proportionate mixing) and a synthetic, community-structured contact network derived from real-world empirical data, loaded as an edgelist |
| **Part 3** | A behavioural response (SIHR with hospitalization-driven contact restriction) and contact-rate heterogeneity (front-line workers vs general public)                                            |

Each notebook writes its outputs to its own subdirectory (`figures/part1/`, `figures/part2/`, `figures/part3/`, and similarly for `data/`), so the three can be re-run independently without clobbering each other's results.

## Installation instructions

1. **Visual Studio Code**

   - Install the latest version of [Visual Studio Code](https://code.visualstudio.com/)
2. **Julia Language**

   - Install the latest version of [Julia](https://julialang.org/downloads/)
   - On Windows, you can also install it through the Microsoft Store
   - Alternatively, use [JuliaUp](https://github.com/JuliaLang/juliaup) for better version management
3. **VS Code Extensions**

   - Open VS Code and install the Julia extension
   - You can find it by searching for "Julia" in the Extensions tab
   - If the Julia extension crashes, you might have to set the path to the Julia executable manually
   - Go to Preferences → Settings → search for "Julia executable path" and paste the path to your `julia.exe`
4. **Verify Installation**

   - Open the Julia REPL in VS Code by pressing:
     - Windows: `Alt+J, Alt+O`
     - macOS: `Ctrl+J, Ctrl+O`
     - Linux: `Ctrl+J, Ctrl+O`
   - Alternatively, use the Command Palette (F1 or Ctrl/Cmd+Shift+P) and search for "Julia: Start REPL"
   - If the REPL opens successfully, your setup is working correctly
5. **Run the Setup**

   - Clone or download this repository
     - Open the source control tab in VS Code and click "Clone repository" (you need git installed for this)
     - Paste the link from Slack and select the path where you want to clone it
     - After successful cloning, you should see all the necessary files in your workspace
   - Open the notebook files (`Part1.ipynb`, `Part2.ipynb`, and `Part3.ipynb`) in VS Code
   - Select your version of Julia as the kernel
   - You can run the first two cells to test if everything works and to already install all necessary packages
   - Allow network access for Julia if prompted

If you encounter any issues during installation, please check the [Julia documentation](https://docs.julialang.org/) or open an issue in this repository. Reinstalling the Julia Extension as well as restarting VS Code (preferably running it as an administrator) might help.

---

## Part 1: Network Types and Epidemic Modelling

Three synthetic networks (random, small-world, preferential attachment) on a homogeneous population. The notebook walks through a single model in detail, then sweeps multiple networks and seeds for statistical comparison.

> **Heads-up:** the *Model Setup* section initializes a single `model` object that is consumed by the next four sections (Network Analysis → Network Visualization → SIR Dynamics → Epidemic Trajectory Visualization). From *Single Run Visualization* onward, each helper function builds its own fresh model. The notebook flags both boundaries explicitly.

### Task 1: Getting Started with Random Networks

1. Start with graph type `:random`.
2. Run the notebook; these steps will run automatically:
   - Generate the network (done by `initialize`)
   - Calculate and visualize the degree distribution
   - Calculate network properties (clustering coefficient, component sizes, diameter)
   - Calculate centrality measures for nodes
   - Run the SIR model on the network with default parameters

### Task 2: Comparing Network Types

1. Repeat the same analysis with graph types `:smallworld` and `:preferential`.
2. Compare degree distributions and other graph-based measures.
3. Compare the dynamics of the outbreaks.
4. Compare the final sizes.
5. Investigate whether there is a relationship between any of the graph measures and the final epidemic size (number of infected individuals at the end of the outbreak).

### Task 3: Parameter Variation

1. Choose one of the three network types.
2. Vary the `mean_degree` and observe the effect on epidemic dynamics and final state.
3. Vary the index case (`patient_zero`) — start the outbreak from a random individual or one with high `degree_centrality`, `betweenness_centrality`, or `eigenvector_centrality` — and observe how it affects the epidemic.

---

## Part 2: Configuration Networks and Real-World Data

Three new network types built from realistic contact data:

- **Configuration model** — built from an empirical degree distribution (POLYMOD-like contact-survey data, available for several countries in `degs/deg_dist_<COUNTRY>.csv`).
- **Proportionate mixing** — degrees sampled from a negative binomial distribution fitted to the same data.
- **Edgelist** — a **synthetic** contact network (N = 1000, mean degree 10) **derived from real-world empirical data** on individuals' social ties (network size, family/friend mix, density), loaded from `degs/edgelist_n1000.csv`. Built by generating families first and then connecting them through friendships, which produces strong community structure and high clustering.

### Task 1: Data Preparation and Analysis

1. Familiarize yourself with the data files in `degs/`.
2. Choose one of the country-specific degree distributions, named `deg_dist_<COUNTRY>.csv`.
3. Calculate and visualize basic statistics of the degree distribution.
4. Fit a negative binomial distribution to the data using maximum likelihood estimation. The output parameters are the mean μ and the dispersion θ.
5. Plot a histogram of the distribution and overlay the fitted PDF; visualize how well the fit captures the observed data.

### Task 2: Configuration Network and Proportionate Mixing

1. Use the file `deg_dist_<COUNTRY>.csv` to generate a configuration network and run the SIR model on it.
2. Compare the results to those from Part 1. How does the structure of this network compare to random / small-world / preferential attachment?
3. **Interpretation:** in which aspects is the generated configuration network *not* a realistic model of the real contact network measured by POLYMOD?
4. Run the proportionate mixing network with the negative-binomial parameters estimated from the data. Are the results similar to the configuration network?
5. Vary the mean degree (and optionally the dispersion) of the negative binomial distribution and observe how this influences network structure and epidemic dynamics.

### Task 3: A Community-Structured Edgelist Network

1. Load the edgelist network from `degs/edgelist_n1000.csv` (this happens automatically near the top of the notebook). Inspect its size, mean degree, clustering coefficient, and community structure.
2. Run the SIR model on this network and compare its trajectory to the configuration and proportionate-mixing networks.
3. Use the comparison plots at the bottom of the notebook (network-metrics bar chart and centrality-distribution boxplots) to characterise how the community-structured edgelist network differs from the degree-distribution-based configuration and proportionate-mixing networks.

---

## Part 3: Behavioural Response and Contact Heterogeneity

A closed community of 1,000 individuals on a small-world network, with two new ingredients:

1. **SIHR dynamics** — infected agents may develop severe disease and be hospitalised (removed from transmission) before recovering.
2. **A behavioural feedback loop** — as hospitalisations rise, individuals adopt contact-restricting behaviour (logistic adoption curve, absorbing state, halves their contacts).
3. **Contact-rate heterogeneity** — a minority of front-line workers (`contact_rate = :high`) drives spread; the rest (`contact_rate = :low`) transmits at `low_risk_factor` of the front-line rate. This dimension is dormant in Sections A and B (where `low_risk_factor = 1.0` makes the two groups equivalent) and switches on in Sections C and D.

> **Heads-up:** Section B is fully tunable from a single source of truth — `behavior_param_sets` near the top of the section. Edit α, β, `max_prob`, and `contact_reduction` there; labels, simulations, and all downstream plots follow automatically.

### Task 1: SIHR Baseline (Section A)

1. Run the plain SIHR model on a small-world network (`mean_degree = 10`, no behavioural response).
2. Compare its trajectory to the SIR version (`use_hospitalization=false`). Where does the hospitalized compartment sit relative to peak infected? Does it materially change the final size?

### Task 2: Hospitalisation-Driven Behavioural Response (Section B)

1. Visualise the logistic adoption curve `p_adopt(s) = max_prob / (1 + e^(-α (s - β)))` for the three default parameter sets.
2. Run the full model (SIHR + behaviour) for the baseline (no behaviour) and the three behaviour sets. Compare mean infected and hospitalised trajectories across seeds.
3. Compare peak hospitalised and final epidemic size across scenarios.
4. **Tune:** edit `behavior_param_sets` (α, β, `max_prob`) and `contact_reduction` to see how earlier / sharper / stronger responses change the outcome.

### Task 3: Contact Heterogeneity — Where Do the High-Contact Workers Sit? (Section C)

1. Switch to a preferential-attachment network and split the population: 10% high-contact front-line workers, 90% lower-contact general public (`low_risk_factor = 0.3`).
2. Compare three placement strategies for the high-contact minority (`high_contact`): `:random`, `:maxdegree`, `:maxbetweenness`.
3. Look at the resulting boxplots (epidemic duration, peak infected, susceptible fraction remaining). Does placing high-contact individuals at hubs really act as a super-spreader effect?

### Task 4: Sensitivity Sweep (Section D)

1. On the same preferential-attachment network (high-contact workers placed at hubs), sweep over a grid of `fraction_high_contact` × `low_risk_factor`.
2. Examine the heatmap of mean peak hospitalised — your proxy for healthcare-system burden — and identify which corner of the parameter space is most dangerous.
3. Discuss: what's more leverage-rich for intervention, *reducing the number* of high-contact individuals or *reducing the transmission* of the lower-contact majority?
