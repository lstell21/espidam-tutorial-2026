"""
    initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, high_contact,
               fraction_high_contact, trans_prob, days_to_recovered, seed, r̂, p̂, low_risk_factor,
               use_hospitalization, hospitalization_prob, days_to_hospital_recovery,
               use_ppe_adoption, max_ppe_reduction, ppe_half_saturation, custom_graph,
               edgelist_path, degrees)

Initialize the model with specified parameters.

# Arguments
- `network_type`: The type of network to create. Can be `:random`, `:smallworld`, `:preferential`, `:configuration`, `:proportionatemixing`, `:edgelist`, or `:custom`.
- `mean_degree`: The mean degree of the network. Default is 4.
- `n_nodes`: The number of nodes in the network. Default is 1000.
- `dispersion`: The dispersion parameter for negative binomial distributions. Default is 0.1.
- `patient_zero`: How to select the initial infected agent. Can be `:random`, `:maxdegree`, `:maxbetweenness`, or `:maxeigenvector`. Default is `:random`.
- `high_contact`: How high-contact-rate agents are distributed. Can be `:random`, `:maxdegree`, `:maxbetweenness`, or `:maxeigenvector`. Default is `:random`.
- `fraction_high_contact`: The fraction of agents with high contact rate (e.g., HCWs). Default is 0.1.
- `trans_prob`: The transmission probability for high-contact agents. Default is 0.1.
- `days_to_recovered`: Days until recovery from infection. Default is 14.
- `seed`: Random seed. Default is 42.
- `r̂`: Negative binomial r parameter (for proportionate mixing). Default is nothing.
- `p̂`: Negative binomial p parameter (for proportionate mixing). Default is nothing.
- `low_risk_factor`: Transmission multiplier for low-contact agents (0–1). Default is 1.0.
- `use_hospitalization`: Whether to enable hospitalization dynamics. Set to `false` for SIR, `true` (default) for SIHR.
- `hospitalization_prob`: Probability of hospitalization for an infected agent (rolled once, on `hospitalization_day`). Default is 0.1.
- `hospitalization_day`: Day of infection on which the hospitalization roll is made (1 = newly infected agents may be hospitalized before their first transmission step). Default is 1.
- `days_to_hospital_recovery`: Days until recovery from hospitalization. Default is 7.
- `use_ppe_adoption`: Whether to enable behavioral PPE-adoption dynamics (HCW transmission is scaled down as the hospitalization burden rises). Default is false.
- `max_ppe_reduction`: Maximum fractional reduction in HCW transmission achievable through PPE (0–1). E.g. 0.8 means PPE can cut HCW transmission by up to 80% at high hospitalization burden. Default is 0.8.
- `ppe_half_saturation`: Hospitalized fraction at which the behavioral response reaches half of `max_ppe_reduction`. Smaller values mean HCWs react to the outbreak earlier. Default is 0.02.
- `custom_graph`: A pre-loaded graph object to use instead of creating a new one. If provided, `n_nodes`/`mean_degree` are inferred from the graph; `network_type` is kept as specified by the caller. Default is nothing.
- `edgelist_path`: Path to the edgelist file when `network_type` is `:edgelist`. Default is "degs/network".
- `degrees`: Degree sequence vector for `:configuration` network type. Default is nothing.

# Returns
- `model`: The initialized agent-based model.
"""
function initialize(;
    network_type::Symbol,
    mean_degree::Integer=4,
    n_nodes::Integer=1000,
    dispersion::Float64=0.1,
    patient_zero::Symbol=:random,
    high_contact::Symbol=:random,
    fraction_high_contact::Float64=0.1,
    trans_prob::Float64=0.1,
    days_to_recovered::Integer=14,
    seed=42,
    r̂=nothing,
    p̂=nothing,
    low_risk_factor::Float64=1.0,
    use_hospitalization::Bool=true,
    hospitalization_prob::Float64=0.1,
    hospitalization_day::Integer=1,
    days_to_hospital_recovery::Integer=7,
    use_ppe_adoption::Bool=false,
    max_ppe_reduction::Float64=0.8,
    ppe_half_saturation::Float64=0.02,
    custom_graph=nothing,
    edgelist_path::String="degs/network",
    degrees=nothing
)
    if !(0 <= low_risk_factor <= 1)
        error("low_risk_factor must be between 0 and 1, got $low_risk_factor")
    end

    rng = Xoshiro(seed)

    # Create or use provided graph
    if custom_graph !== nothing
        graph = custom_graph
        n_nodes = nv(graph)
        mean_degree = round(Int, 2 * ne(graph) / nv(graph))
    else
        graph = create_graph(; network_type, mean_degree, n_nodes, dispersion, r̂, p̂, edgelist_path, degrees, rng)
        if network_type == :edgelist
            n_nodes = nv(graph)
            mean_degree = round(Int, 2 * ne(graph) / nv(graph))
        end
    end

    space = GraphSpace(graph)
    properties = create_properties(graph, network_type, n_nodes, mean_degree, dispersion,
                                   patient_zero, high_contact, fraction_high_contact, trans_prob,
                                   days_to_recovered, low_risk_factor, r̂, p̂,
                                   use_hospitalization, hospitalization_prob, hospitalization_day,
                                   days_to_hospital_recovery,
                                   use_ppe_adoption, max_ppe_reduction, ppe_half_saturation)
    model = StandardABM(Person, space; agent_step!, model_step!, properties, rng)
    populate(model, high_contact, fraction_high_contact)
    set_patient_zero!(model, patient_zero)
    return model
end



############################### helper functions ###############################

function create_properties(graph, network_type, n_nodes, mean_degree, dispersion, patient_zero,
                           high_contact, fraction_high_contact, trans_prob, days_to_recovered,
                           low_risk_factor=1.0, r̂=nothing, p̂=nothing,
                           use_hospitalization=true, hospitalization_prob=0.1,
                           hospitalization_day=1, days_to_hospital_recovery=7,
                           use_ppe_adoption=false, max_ppe_reduction=0.8, ppe_half_saturation=0.02)
    low_risk_factor = clamp(low_risk_factor, 0.0, 1.0)

    properties = Dict(
        :graph => graph,
        :network_type => network_type,
        :n_nodes => n_nodes,
        :mean_degree => mean_degree,
        :dispersion => dispersion,
        :patient_zero => patient_zero,
        :high_contact => high_contact,
        :fraction_high_contact => fraction_high_contact,
        :trans_prob => trans_prob,
        :days_to_recovered => days_to_recovered,
        :low_risk_factor => low_risk_factor,
        :use_hospitalization => use_hospitalization,
        :hospitalization_prob => hospitalization_prob,
        :hospitalization_day => hospitalization_day,
        :days_to_hospital_recovery => days_to_hospital_recovery,
        :use_ppe_adoption => use_ppe_adoption,
        :max_ppe_reduction => max_ppe_reduction,
        :ppe_half_saturation => ppe_half_saturation,
        :susceptible_count => n_nodes,
        :infected_count => 1,
        :hospitalized_count => 0,
        :recovered_count => 0)

    if r̂ !== nothing
        properties[:r̂] = r̂
    end
    if p̂ !== nothing
        properties[:p̂] = p̂
    end

    return properties
end

"""
populate(model, high_contact::Symbol, fraction_high_contact::Float64)

Populates the model with agents based on the specified contact-rate distribution.

## Arguments
- `model::AgentBasedModel`: The model to populate with agents.
- `high_contact::Symbol`: How to assign high-contact-rate agents. Possible values are `:random`, `:maxdegree`, `:maxbetweenness`, and `:maxeigenvector`.
- `fraction_high_contact::Float64`: The fraction of agents to assign as high-contact-rate. Default is 0.1.
"""
function populate(model::AgentBasedModel, high_contact::Symbol, fraction_high_contact::Float64=0.1)
    if high_contact == :random
        n_high = Int(round(fraction_high_contact * model.n_nodes))
        n_low  = model.n_nodes - n_high

        for _ in 1:n_high
            add_agent_single!(model, :S, 0, :high)
        end
        for _ in 1:n_low
            add_agent_single!(model, :S, 0, :low)
        end
    else
        centrality_func = if high_contact == :maxdegree
            degree_centrality
        elseif high_contact == :maxbetweenness
            betweenness_centrality
        elseif high_contact == :maxeigenvector
            eigenvector_centrality
        else
            error("Unknown high_contact strategy: $high_contact")
        end

        sorted_nodes = sortperm(centrality_func(model.graph), rev=true)
        n_high = Int(floor(fraction_high_contact * length(sorted_nodes)))
        selected_positions = sorted_nodes[1:n_high]

        # Place high-contact agents at their hub positions FIRST so add_agent_single!
        # for low-contact agents can only land on the remaining empty nodes.
        for pos in selected_positions
            add_agent!(pos, model, :S, 0, :high)
        end
        for _ in 1:(model.n_nodes - n_high)
            add_agent_single!(model, :S, 0, :low)
        end
    end
end

"""
set_patient_zero!(model::AgentBasedModel, patient_zero::Symbol)

Set the initial infected agent in the model based on the given `patient_zero` strategy.

## Arguments
- `model::AgentBasedModel`: The agent-based model.
- `patient_zero::Symbol`: The strategy to determine the initial infected agent. Possible values are `:random`, `:maxdegree`, `:maxbetweenness`, and `:maxeigenvector`.
"""
function set_patient_zero!(model::AgentBasedModel, patient_zero::Symbol)
    if patient_zero == :random
        random_agent(model).status = :I
        return
    end
    node = if patient_zero == :maxdegree
        argmax(degree_centrality(model.graph))
    elseif patient_zero == :maxbetweenness
        argmax(betweenness_centrality(model.graph))
    elseif patient_zero == :maxeigenvector
        argmax(eigenvector_centrality(model.graph))
    else
        error("Unknown patient_zero strategy: $patient_zero")
    end
    first(agents_in_position(node, model)).status = :I
end
