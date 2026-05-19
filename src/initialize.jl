"""
    initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_risk,
               fraction_high_risk, trans_prob, days_to_recovered, seed, r̂, p̂, low_risk_factor,
               use_hospitalization, hospitalization_prob, days_to_hospital_recovery,
               custom_graph, edgelist_path)

Initialize the model with specified parameters.

# Arguments
- `network_type`: The type of network to create. Can be `:random`, `:smallworld`, `:preferential`, `:configuration`, `:proportionatemixing`, `:edgelist`, or `:custom`.
- `mean_degree`: The mean degree of the network. Default is 4.
- `n_nodes`: The number of nodes in the network. Default is 1000.
- `dispersion`: The dispersion parameter for negative binomial distributions. Default is 0.1.
- `patient_zero`: How to select the initial infected agent. Can be `:random`, `:maxdegree`, `:maxbetweenness`, or `:maxeigenvector`. Default is `:random`.
- `high_risk`: How high-risk agents are distributed. Can be `:random`, `:maxdegree`, `:maxbetweenness`, or `:maxeigenvector`. Default is `:random`.
- `fraction_high_risk`: The fraction of high-risk agents. Default is 0.1.
- `trans_prob`: The transmission probability. Default is 0.1.
- `days_to_recovered`: Days until recovery from infection. Default is 14.
- `seed`: Random seed. Default is 42.
- `r̂`: Negative binomial r parameter (for proportionate mixing). Default is nothing.
- `p̂`: Negative binomial p parameter (for proportionate mixing). Default is nothing.
- `low_risk_factor`: Transmission multiplier for low-risk agents. Default is 1.0.
- `use_hospitalization`: Whether to enable hospitalization dynamics. Set to `false` for SIR (Part 1), `true` (default) for SIHR (Part 2).
- `hospitalization_prob`: Probability of hospitalization for an infected agent. Default is 0.1.
- `days_to_hospital_recovery`: Days until recovery from hospitalization. Default is 7.
- `custom_graph`: A pre-loaded graph object to use instead of creating a new one. If provided, `network_type` is automatically set to `:custom` and `n_nodes`/`mean_degree` are inferred. Default is nothing.
- `edgelist_path`: Path to the edgelist file when `network_type` is `:edgelist`. Default is "degs/network".

# Returns
- `model`: The initialized agent-based model.
"""
function initialize(;
    network_type::Symbol,
    mean_degree::Integer=4,
    n_nodes::Integer=1000,
    dispersion::Float64=0.1,
    patient_zero::Symbol=:random,
    high_risk::Symbol=:random,
    fraction_high_risk::Float64=0.1,
    trans_prob::Float64=0.1,
    days_to_recovered::Integer=14,
    seed=42,
    r̂=nothing,
    p̂=nothing,
    low_risk_factor::Float64=1.0,
    use_hospitalization::Bool=true,
    hospitalization_prob::Float64=0.1,
    days_to_hospital_recovery::Integer=7,
    custom_graph=nothing,
    edgelist_path::String="degs/network"
)
    # Validate low_risk_factor
    if !(0 <= low_risk_factor <= 1)
        error("low_risk_factor must be between 0 and 1, got $low_risk_factor")
    end

    # Create or use provided graph
    if custom_graph !== nothing
        # Use the provided custom graph
        graph = custom_graph
        network_type = :custom
        n_nodes = nv(graph)
        mean_degree = round(Int, 2 * ne(graph) / nv(graph))
    else
        # Create a new graph
        graph = create_graph(; network_type, mean_degree, n_nodes, dispersion, r̂, p̂, edgelist_path)
        # Update n_nodes and mean_degree for edgelist graphs
        if network_type == :edgelist
            n_nodes = nv(graph)
            mean_degree = round(Int, 2 * ne(graph) / nv(graph))
        end
    end

    space = GraphSpace(graph)
    # set up properties
    properties = create_properties(graph, network_type, n_nodes, mean_degree, dispersion,
                                   patient_zero, high_risk, fraction_high_risk, trans_prob,
                                   days_to_recovered, low_risk_factor, r̂, p̂,
                                   use_hospitalization, hospitalization_prob, days_to_hospital_recovery)
    # set up RNG
    rng = Xoshiro(seed)
    # create the model
    model = StandardABM(Person, space; agent_step!, model_step!, properties, rng)
    # add agents, if high_risk is random, add high risk agents randomly
    populate(model, high_risk, fraction_high_risk)
    set_patient_zero!(model, patient_zero)
    return model
end



############################### helper functions ###############################

"""
create_properties(graph, network_type, n_nodes, mean_degree, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, days_to_recovered, low_risk_factor, r̂, p̂)

Create a dictionary of properties for the simulation.

# Arguments
- `graph`: The graph object representing the network structure.
- `network_type`: The type of network.
- `n_nodes`: The number of nodes in the network.
- `mean_degree`: The mean degree of the network.
- `dispersion`: The dispersion parameter for the network.
- `patient_zero`: The initial infected node.
- `high_risk`: Whether high-risk nodes are present in the network.
- `fraction_high_risk`: The fraction of high-risk nodes in the network.
- `trans_prob`: The transmission probability.
- `days_to_recovered`: The number of days it takes for an infected node to recover.
- `low_risk_factor`: Factor to multiply the transmission probability for low risk agents.
- `r̂`: The r parameter for negative binomial distribution.
- `p̂`: The p parameter for negative binomial distribution.

# Returns
- `properties`: A dictionary containing the properties for the simulation.

"""
function create_properties(graph, network_type, n_nodes, mean_degree, dispersion, patient_zero, high_risk, fraction_high_risk, trans_prob, days_to_recovered, low_risk_factor=1.0, r̂=nothing, p̂=nothing, use_hospitalization=true, hospitalization_prob=0.1, days_to_hospital_recovery=7)
    # Ensure low_risk_factor is between 0 and 1
    low_risk_factor = clamp(low_risk_factor, 0.0, 1.0)
    
    properties = Dict(
        :graph => graph,
        :network_type => network_type,
        :n_nodes => n_nodes,
        :mean_degree => mean_degree,
        :dispersion => dispersion,
        :patient_zero => patient_zero,
        :high_risk => high_risk,
        :fraction_high_risk => fraction_high_risk,
        :trans_prob => trans_prob,
        :days_to_recovered => days_to_recovered,
        :low_risk_factor => low_risk_factor,
        :use_hospitalization => use_hospitalization,
        :hospitalization_prob => hospitalization_prob,
        :days_to_hospital_recovery => days_to_hospital_recovery,
        :susceptible_count => n_nodes,
        :infected_count => 1,
        :hospitalized_count => 0,
        :recovered_count => 0)
    
    # Add r̂ and p̂ to properties if provided
    if r̂ !== nothing
        properties[:r̂] = r̂
    end
    if p̂ !== nothing
        properties[:p̂] = p̂
    end
    
    return properties
end

"""
populate(model, high_risk::Symbol, fraction_high_risk::Float64)

Populates the model with agents based on the specified risk distribution.

## Arguments
- `model::AgentBasedModel`: The model to populate with agents.
- `high_risk::Symbol`: The risk distribution to use. Possible values are `:random`, `:maxdegree`, `:maxbetweenness`, and `:maxeigenvector`.
- `fraction_high_risk::Float64`: The fraction of nodes to assign as high-risk agents. Default is 0.1.

## Details
- If `high_risk` is `:random`, agents are randomly assigned as high-risk or low-risk.
- If `high_risk` is `:maxdegree`, fraction_high_risk agents with the highest degree centrality are assigned as high-risk.
- If `high_risk` is `:maxbetweenness`, fraction_high_risk agents with the highest betweenness centrality are assigned as high-risk.
- If `high_risk` is `:maxeigenvector`, fraction_high_risk agents with the highest eigenvector centrality are assigned as high-risk.
"""
function populate(model::AgentBasedModel, high_risk::Symbol, fraction_high_risk::Float64=0.1)
    if high_risk == :random
        n_high_risk = Int(round(fraction_high_risk * model.n_nodes))
        n_low_risk = model.n_nodes - n_high_risk

        for _ in 1:n_high_risk
            add_agent_single!(model, :S, 0, :high)
        end
        for _ in 1:n_low_risk
            add_agent_single!(model, :S, 0, :low)
        end
    else
        centrality_func = if high_risk == :maxdegree
            degree_centrality
        elseif high_risk == :maxbetweenness
            betweenness_centrality
        elseif high_risk == :maxeigenvector
            eigenvector_centrality
        else
            error("Unknown high_risk strategy: $high_risk")
        end

        sorted_nodes = sortperm(centrality_func(model.graph), rev=true)
        selected_positions = sorted_nodes[1:Int(floor(fraction_high_risk * length(sorted_nodes)))]
        for i in 1:model.n_nodes
            if i in selected_positions
                add_agent!(i, model, :S, 0, :high)
            else
                add_agent_single!(model, :S, 0, :low)
            end
        end
    end
end

"""
set_patient_zero!(model::AgentBasedModel, patient_zero::Symbol)

Set the initial infected agent in the model based on the given `patient_zero` strategy.

## Arguments
- `model::AgentBasedModel`: The agent-based model.
- `patient_zero::Symbol`: The strategy to determine the initial infected agent. Possible values are `:random`, `:maxdegree`, `:maxbetweenness`, and `:maxeigenvector`.

## Details
- If `patient_zero` is `:random`, a random agent in the model will be set as infected.
- If `patient_zero` is `:maxdegree`, the agent with the highest degree centrality in the model's graph will be set as infected.
- If `patient_zero` is `:maxbetweenness`, the agent with the highest betweenness centrality in the model's graph will be set as infected.
- If `patient_zero` is `:maxeigenvector`, the agent with the highest eigenvector centrality in the model's graph will be set as infected.
"""
function set_patient_zero!(model::AgentBasedModel, patient_zero::Symbol)
    if patient_zero == :random
        random_agent(model).status = :I
    elseif patient_zero == :maxdegree
        model[argmax(degree_centrality(model.graph))].status = :I
    elseif patient_zero == :maxbetweenness
        model[argmax(betweenness_centrality(model.graph))].status = :I
    elseif patient_zero == :maxeigenvector
        model[argmax(eigenvector_centrality(model.graph))].status = :I
    end
end