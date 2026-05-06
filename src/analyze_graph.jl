"""
    analyze_graph(g::AbstractGraph)

Analyze a graph and compute various network metrics and properties.

This function calculates a comprehensive set of graph metrics including:
- Basic statistics (mean degree, density, clustering coefficient)
- Centrality measures (degree, betweenness, closeness, eigenvector)
- Component analysis
- Clique information

# Arguments
- `g::AbstractGraph`: The graph to analyze

# Returns
A dictionary containing the following keys:
- `"summary"`: DataFrame with scalar network metrics
- `"centrality"`: DataFrame with node-level centrality measures
- `"degree_distribution"`: Dictionary mapping degrees to their frequency
- `"connected_components"`: Vector of components (each a vector of vertex IDs)
- `"component_lengths"`: Vector of component sizes
- `"maximal_cliques"`: Vector of maximal cliques in the graph

# Examples
```julia
using Graphs, DataFrames
g = erdos_renyi(100, 0.1)
results = analyze_graph(g)
summary_metrics = results["summary"]
centrality_data = results["centrality"]
```
"""
function analyze_graph(g::AbstractGraph)
    # Calculate network metrics
    mean_degree = mean(Graphs.degree(g))
    density = Graphs.density(g)
    clustering_coefficient = global_clustering_coefficient(g)
    assortativity = Graphs.assortativity(g)
    diam = is_connected(g) ? diameter(g) : "Graph is not connected"
    degree_distribution = degree_histogram(g)
    
    # Calculate centrality measures
    dg_c = degree_centrality(g)
    btwn_c = betweenness_centrality(g)
    clns_c = closeness_centrality(g)
    eig_c = eigenvector_centrality(g)
    
    # Calculate component information
    cnct_components = connected_components(g)
    comp_lengths = map(length, cnct_components)
    max_comp_length = maximum(comp_lengths)
    max_cliques = maximal_cliques(g)
    
    # Create summary dataframe with scalar values
    results_summary = DataFrame(
        metric = [
            "Mean Degree", 
            "Density", 
            "Clustering Coefficient",
            "Assortativity",
            "Diameter",
            "Number of Connected Components",
            "Max Component Length"
        ],
        value = [
            mean_degree, 
            density, 
            clustering_coefficient,
            assortativity,
            diam,
            length(cnct_components),
            max_comp_length
        ]
    )
    
    # Store more detailed node-level metrics in a separate dataframe
    centrality_measures = DataFrame(
        node_id = 1:nv(g),
        degree_centrality = dg_c,
        betweenness_centrality = btwn_c,
        closeness_centrality = clns_c,
        eigenvector_centrality = eig_c
    )
    
    # Create a dictionary to hold all results
    full_results = Dict(
        "summary" => results_summary,
        "centrality" => centrality_measures,
        "degree_distribution" => degree_distribution,
        "connected_components" => cnct_components,
        "component_lengths" => comp_lengths,
        "maximal_cliques" => max_cliques
    )
    
    return full_results
end

"""
    print_graph_analysis(analysis_results)

Display the graph analysis results in a notebook-friendly format.

# Arguments
- `analysis_results`: The dictionary returned by `analyze_graph`.
"""
function print_graph_analysis(analysis_results)
    # Create a result that will display nicely in Jupyter
    summary_df = analysis_results["summary"]
    
    # Get centrality statistics for display
    centrality_summary = describe(analysis_results["centrality"][:, 2:end])
    
    # Get degree distribution info
    deg_dist = analysis_results["degree_distribution"]
    deg_keys = collect(keys(deg_dist))
    deg_values = collect(values(deg_dist))
    
    deg_info = DataFrame(
        metric = ["Min Degree", "Max Degree", "Most Common Degree"],
        value = [
            minimum(deg_keys),
            maximum(deg_keys),
            deg_keys[findmax(deg_values)[2]]
        ]
    )
    
    # Create a component sizes table and summary
    comp_sizes = analysis_results["component_lengths"]
    component_sizes_df = DataFrame(component_id = 1:length(comp_sizes), size = comp_sizes)
    
    # Create a component size summary table
    size_counts = countmap(comp_sizes)
    component_summary = DataFrame(
        component_size = collect(keys(size_counts)),
        count = collect(values(size_counts))
    )
    sort!(component_summary, :component_size)
    
    # Return a named tuple of dataframes that will display well in Jupyter
    return (
        summary = summary_df,
        centrality = centrality_summary,
        degree_info = deg_info,
        component_sizes = component_sizes_df,
        component_summary = component_summary
    )
end