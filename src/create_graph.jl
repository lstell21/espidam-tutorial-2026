using Random
using CSV, DataFrames
"""
    create_graph(; network_type, mean_degree, n_nodes, dispersion = 0.1, β = 0.1, k = 3, r̂ = nothing, p̂ = nothing, edgelist_path = nothing, degrees = nothing)

Create a graph based on the specified network type.

# Arguments
- `network_type`: The type of network to create. Can be `:random`, `:smallworld`, `:preferential`, `:configuration`, `:proportionate`, or `:edgelist`.
- `n_nodes`: The number of nodes in the graph. For :configuration, the number of nodes is fixed to 1000. Default is 1000.
- `mean_degree`: The average degree of the nodes in the graph. For :preferential, k is used as mean_degree/2 (for even numbers). Default is 4.
- `dispersion`: The dispersion parameter for the negative binomial distribution, used only when `network_type` is `:proportionate` and r̂ and p̂ are not provided. Default is 0.1.
- `β`: The rewiring probability, used only when `network_type` is `:smallworld`. Default is 0.1.
- `k`: The number of edges to attach from a new node to existing nodes, used only when `network_type` is `:preferential`. Default is 3.
- `r̂`: The r parameter for negative binomial distribution, used only when `network_type` is `:proportionate`. If not provided, will be calculated from mean_degree and dispersion.
- `p̂`: The p parameter for negative binomial distribution, used only when `network_type` is `:proportionate`. If not provided, will be calculated from mean_degree and dispersion.
- `edgelist_path`: The path to the edgelist file, used only when `network_type` is `:edgelist`. The file is a CSV with one edge per row and two node-id columns (a header row is optional). Default is "degs/edgelist_n1000.csv".
- `degrees`: Degree sequence vector (required for `:configuration` network type).

# Returns
- `graph`: The created graph.

# Examples
```julia
g = create_graph(; network_type = :random,  mean_degree = 4)
g = create_graph(; network_type = :proportionate, n_nodes = 1000, r̂ = 5.0, p̂ = 0.4)
g = create_graph(; network_type = :edgelist, edgelist_path = "degs/edgelist_n1000.csv")
g = create_graph(; network_type = :configuration, degrees = my_degree_vector)
```
"""
function create_graph(; network_type::Symbol, mean_degree::Integer, n_nodes::Integer=1000, dispersion::Float64=0.1, β::Float64=0.1, k::Integer=4, r̂=nothing, p̂=nothing, edgelist_path::String="degs/edgelist_n1000.csv", degrees=nothing, rng::AbstractRNG=Random.default_rng())
    if network_type == :random
        graph = Graphs.erdos_renyi(n_nodes, mean_degree / n_nodes; rng)
    elseif network_type == :smallworld
        graph = Graphs.newman_watts_strogatz(n_nodes, mean_degree, β::Float64; rng)
    elseif network_type == :preferential
        graph = Graphs.barabasi_albert(n_nodes, Int(round(mean_degree / 2)); rng)
    elseif network_type == :configuration
        if degrees === nothing
            error("For network_type=:configuration, the `degrees` argument must be provided (a vector of node degrees).")
        end
        graph = Graphs.random_configuration_model(1000, degrees; rng)
    elseif network_type == :proportionatemixing || network_type == :proportionate
        if isnothing(r̂) || isnothing(p̂)
            error("Both r̂ and p̂ must be provided for proportionate mixing networks")
        end

        # Generate degree sequence from negative binomial distribution
        degree_sequence = rand(rng, NegativeBinomial(r̂, p̂), n_nodes)
        
        # Ensure all degrees are within valid bounds [0, n_nodes-1]
        degree_sequence = clamp.(degree_sequence, 0, n_nodes - 1)
        
        # No graphicality or even-sum check is needed here: expected_degree_graph (the
        # Chung–Lu model) treats the sequence as *expected* degrees and wires each pair of
        # nodes independently, so any non-negative sequence is admissible. (Those checks
        # belong to random_configuration_model, used for the :configuration network.)
        # Build the proportionate-mixing graph directly from the NB degree sequence,
        # threading the model rng so wiring is reproducible like every other generator here.
        graph = expected_degree_graph(degree_sequence; rng)
    elseif network_type == :edgelist
        # Load an undirected contact network from a CSV edgelist: one edge per row,
        # the first two columns are the node ids of each edge's endpoints.
        edge_df = CSV.read(edgelist_path, DataFrame)
        n = max(maximum(edge_df[!, 1]), maximum(edge_df[!, 2]))
        graph = SimpleGraph(n)
        for row in eachrow(edge_df)
            add_edge!(graph, row[1], row[2])
        end
    elseif network_type == :custom
        # Custom network type is handled differently in initialize.jl
        # This should not be called directly through create_graph
        error("Custom network type should use the custom_graph parameter in initialize() instead of create_graph()")
    else
        error("Unknown network type: $network_type")
    end
    return graph
end