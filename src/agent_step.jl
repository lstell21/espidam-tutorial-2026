"""
    agent_step!(person::Person, model::AgentBasedModel)

Update the state of an agent in the ABM model for one time step.

# Arguments
- `person::Person`: The agent whose state needs to be updated.
- `model::AgentBasedModel`: The ABM model containing the agent.

# Description
- If infected, increment days infected and recover if enough time has passed
- If infected, potentially transmit to nearby susceptible agents with probability
  adjusted by risk level
- Handle hospitalization transitions for infected agents (if `model.use_hospitalization` is true)
  - When `use_hospitalization=true` (default): Some infected agents may be hospitalized on day 1
  - When `use_hospitalization=false`: No hospitalization occurs (used in Part 1)
"""
function agent_step!(person::Person, model::AgentBasedModel)
    # Check if hospitalization is enabled (defaults to true for backward compatibility)
    use_hospitalization = model.use_hospitalization
    
    # Handle recovery if infected
    if person.status == :I
        person.days_infected += 1
        
        # Check if agent should be hospitalized (only on day 1 of infection)
        if use_hospitalization && person.days_infected == 1 && rand() < model.hospitalization_prob
            person.status = :H
            return
        end
        
        # Check if agent recovers from infection
        if person.days_infected >= model.days_to_recovered
            person.status = :R
            return
        end
        
        # Handle infection of neighbors
        for neighbor in nearby_agents(person, model, 1)
            if neighbor.status == :S
                # Calculate transmission probability based on risk
                trans_prob = neighbor.risk == :high ? model.trans_prob : model.trans_prob * model.low_risk_factor
                if rand(abmrng(model)) < trans_prob
                    neighbor.status = :I
                    neighbor.days_infected = 0
                end
            end
        end
        
    elseif use_hospitalization && person.status == :H
        # Increment days infected (hospitalized agents are still infected)
        person.days_infected += 1
        
        # Check if agent recovers from hospital
        if person.days_infected >= model.days_to_hospital_recovery
            person.status = :R
        end
        
        # Hospitalized agents don't transmit (isolated)
    end
end
