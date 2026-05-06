using GraphIO.EdgeList
"""
    create_graph(; network_type, mean_degree, n_nodes, dispersion = 0.1, β = 0.1, k = 3, r̂ = nothing, p̂ = nothing, edgelist_path = nothing)

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
- `edgelist_path`: The path to the edgelist file, used only when `network_type` is `:edgelist`. Default is "degs/network".

# Returns
- `graph`: The created graph.

# Examples
```julia
g = create_graph(; network_type = :random,  mean_degree = 4)
g = create_graph(; network_type = :proportionate, n_nodes = 1000, r̂ = 5.0, p̂ = 0.4)
g = create_graph(; network_type = :edgelist, edgelist_path = "degs/network")
```
"""
function create_graph(; network_type::Symbol, mean_degree::Integer, n_nodes::Integer=1000, dispersion::Float64=0.1, β::Float64=0.1, k::Integer=4, r̂=nothing, p̂=nothing, edgelist_path::String="degs/network")
    if network_type == :random
        graph = Graphs.erdos_renyi(n_nodes, mean_degree / n_nodes)
    elseif network_type == :smallworld
        graph = Graphs.newman_watts_strogatz(n_nodes, mean_degree, β::Float64)
    elseif network_type == :preferential
        graph = Graphs.barabasi_albert(n_nodes, Int(round(mean_degree / 2)))
    elseif network_type == :configuration
        graph = Graphs.random_configuration_model(1000, degrees[!, 1])
    elseif network_type == :proportionatemixing || network_type == :proportionate
        if isnothing(r̂) || isnothing(p̂)
            error("Both r̂ and p̂ must be provided for proportionate mixing networks")
        end
        
        # Generate degree sequence from negative binomial distribution
        degree_sequence = rand(NegativeBinomial(r̂, p̂), n_nodes)
        
        # Ensure all degrees are within valid bounds [0, n_nodes-1]
        degree_sequence = clamp.(degree_sequence, 0, n_nodes - 1)
        
        # Ensure the sum is even for a valid graph
        if sum(degree_sequence) % 2 != 0
            degree_sequence[1] = max(0, degree_sequence[1] - 1)
        end
        
        # Validate that the degree sequence is graphical
        if !Graphs.isgraphical(degree_sequence)
            # If not graphical, fall back to a configuration model with fixed mean degree
            @warn "Generated degree sequence is not graphical for n_nodes=$n_nodes, falling back to configuration model"
            degree_sequence = fill(mean_degree, n_nodes)
            # Adjust to make sum even
            if sum(degree_sequence) % 2 != 0
                degree_sequence[end] += 1
            end
        end
        
        # Create the graph with the validated degree sequence
        graph = expected_degree_graph(degree_sequence)
    elseif network_type == :edgelist
        # Load graph from edgelist file
        graph = loadgraph(edgelist_path, "graph_key", EdgeListFormat())
        
        # Remove any double edges
        for e in edges(graph)
            if has_edge(graph, src(e), dst(e))
                rem_edge!(graph, src(e), dst(e))
            end
        end
        
        # Convert to undirected simple graph
        graph = SimpleGraph(graph)
    elseif network_type == :custom
        # Custom network type is handled differently in initialize.jl
        # This should not be called directly through create_graph
        error("Custom network type should use the custom_graph parameter in initialize() instead of create_graph()")
    else
        error("Unknown network type: $network_type")
    end
    return graph
end