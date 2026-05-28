"""
    run_simulations(; network_type, mean_degree, patient_zero, high_contact, fraction_high_contact, ...)

Run simulations for an epidemiological model with different combinations of parameters.

# Arguments
- `network_type::Symbol`: The type of network to use for the simulation.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Default is `:random`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. Default is 1000.
- `dispersion::Float64`: The dispersion of the network. Default is 0.1.
- `high_contact::Symbol`: How high-contact-rate agents are distributed. Default is `:random`.
- `fraction_high_contact::Float64`: The fraction of high-contact-rate agents. Default is 1.0.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.
- `r̂`: The r parameter for negative binomial distribution. Default is nothing.
- `p̂`: The p parameter for negative binomial distribution. Default is nothing.
- `low_risk_factor::Float64`: Transmission multiplier for low-contact agents (0–1). Default is 1.0.
- `use_hospitalization::Bool`: Whether to enable hospitalization. Default is true.
- `use_ppe_adoption::Bool`: Whether to enable behavioral PPE adoption dynamics. Default is false.
- `max_ppe_reduction::Float64`: Maximum fractional reduction in HCW transmission via PPE (0–1). Default is 0.8.
- `ppe_half_saturation::Float64`: Hospitalized fraction at which PPE response reaches half its maximum. Default is 0.02.
- `degrees`: Degree sequence vector for `:configuration` network type. Default is nothing.

# Returns
- `mdf::DataFrame`: A DataFrame containing the simulation results.
"""
function run_simulations(; network_type::Symbol, mean_degree::Int, n_nodes::Int=1000,
                        dispersion::Float64=0.1, patient_zero::Symbol=:random,
                        high_contact::Symbol=:random, fraction_high_contact::Float64=1.0,
                        low_risk_factor::Float64=1.0, trans_prob::Float64=0.1,
                        n_steps::Int=100, r̂=nothing, p̂=nothing,
                        use_hospitalization::Bool=true,
                        hospitalization_prob::Float64=0.1,
                        hospitalization_day::Int=1,
                        days_to_hospital_recovery::Int=7,
                        use_ppe_adoption::Bool=false, max_ppe_reduction::Float64=0.8,
                        ppe_half_saturation::Float64=0.02,
                        degrees=nothing)
    if !(0 <= low_risk_factor <= 1)
        error("low_risk_factor must be between 0 and 1, got $low_risk_factor")
    end

    parameters = Dict(
        :seed => rand(UInt16, 100),
        :network_type => network_type,
        :mean_degree => mean_degree,
        :n_nodes => n_nodes,
        :dispersion => dispersion,
        :patient_zero => patient_zero,
        :high_contact => high_contact,
        :trans_prob => trans_prob,
        :fraction_high_contact => fraction_high_contact,
        :low_risk_factor => low_risk_factor,
        :days_to_recovered => 14,
        :use_hospitalization => use_hospitalization,
        :hospitalization_prob => hospitalization_prob,
        :hospitalization_day => hospitalization_day,
        :days_to_hospital_recovery => days_to_hospital_recovery,
        :use_ppe_adoption => use_ppe_adoption,
        :max_ppe_reduction => max_ppe_reduction,
        :ppe_half_saturation => ppe_half_saturation
    )

    if r̂ !== nothing
        parameters[:r̂] = r̂
    end
    if p̂ !== nothing
        parameters[:p̂] = p̂
    end
    if degrees !== nothing
        parameters[:degrees] = [degrees]  # wrap so paramscan treats the whole sequence as one option
    end

    adata = [:status]
    if use_hospitalization
        mdata = [:susceptible_count, :infected_count, :hospitalized_count, :recovered_count]
    else
        mdata = [:susceptible_count, :infected_count, :recovered_count]
    end

    _, mdf = paramscan(
        parameters,
        initialize;
        mdata=mdata,
        n=n_steps,
        showprogress=false
    )
    return mdf
end
