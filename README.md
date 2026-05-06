# ESPIDAM Tutorial 2026

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

5. **Run the Tutorial**
   - Clone or download this repository
     - Open the source control tab in VS Code and click "Clone repository" (you need git installed for this)
     - Paste the link from Slack and select the path where you want to clone it
     - After successful cloning, you should see all the necessary files in your workspace
   - Open the notebook files (Part1.jpynb and Part2.jpynb) in VS Code
   - Select your version of Julia as the kernel
   - You can run the first two cells to test if everything works and to already install all necessary packages
   - Allow network access for Julia if prompted

If you encounter any issues during installation, please check the [Julia documentation](https://docs.julialang.org/) or open an issue in this repository. Reinstalling the Julia Extension as well as restarting VS Code (preferably running it as an administrator) might help.

## Part 1: Network Types and Epidemic Modeling

### Task 1: Getting Started with Random Networks

1. Start with graph type: `random`
2. Run the notebook; these steps will be done automatically:
   - Generate network (done with initialization)
   - Calculate and visualize degree distribution
   - Calculate network properties (clustering coefficient, component sizes, diameter)
   - Calculate centrality measures for nodes
   - Run SIR model on the network with default parameters
3. Start with one graph type to go through all steps

### Task 2: Comparing Network Types

1. Repeat the same analysis with graph types `smallworld` and `preferential`
2. Compare degree distributions and other graph-based measures
3. Compare the dynamics of the outbreaks
4. Compare the final sizes
5. Investigate whether there is a relationship between any of the graph measures and the final epidemic size (number of infected individuals at the end of the outbreak)

### Task 3: Parameter Variation

1. Choose one of the 3 network types
2. Vary the mean degree and observe the effect on epidemic dynamics and final state
3. Vary the index case (start the outbreak from an individual with high degree_centrality, high betweenness_centrality or high eigenvector_centrality) and observe how it affects the epidemic

---

## Part 2: Configuration Networks and Real-World Data

We will explore the network types `configuration network` and `proportionate mixing network`. We do this with a degree distribution based on data, and with a degree distribution sampled from a negative binomial distribution.

### Task 1: Data Preparation and Analysis

1. Make yourself familiar with the data files.
2. Choose one of the country-specific degree distributions, named `deg_dist_"COUNTRY".csv`
3. Calculate and visualize the basic statistics of the degree distribution
4. Fit a negative binomial distribution to the data using maximum likelihood estimation. The output parameters are the mean mu and the dispersion theta.
5. Plot a histogram of the distribution. Visualize how well the fitted distribution captures the observed data

### Task 2: Configuration Network and Proportionate Mixing

1. You now use the file `deg_dist_"COUNTRY".csv` to generate a configuration network and run the model on this network. 
2. Compare the results to those from the previous practical. How does the structure of this network compare to the other network types? 
3. Interpretation: in which aspects is the generated network not a realistic model of the real contact network measured by POLYMOD?
4. Run also the proportionate mixing network with the parameters for the negative binomial distribution estimated (roughly) from the data. Are the results similar to the results from the configuration network?
5. Now vary the mean degree of the negative binomial distribution and observe how this influences network structure and epidemic dynamics. (Optional: vary the dispersion parameter and observe how it influences the dynamics). 
6. How would you include the impact of mask use in this model? One possibility is that the transmission probability is reduced by a certain factor `0<m<1`. Vary the factor m and observe the impact on the epidemic outcome. 
7. Does it matter for the effectiveness of mask use whether the index case is a random individual or has high centrality?