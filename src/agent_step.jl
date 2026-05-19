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
    if person.status == :I
        person.days_infected += 1

        # Check if agent should be hospitalized (only on day 1 of infection)
        if model.use_hospitalization && person.days_infected == 1 && rand() < model.hospitalization_prob
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
                trans_prob = neighbor.risk == :high ? model.trans_prob : model.trans_prob * model.low_risk_factor
                if rand(abmrng(model)) < trans_prob
                    neighbor.status = :I
                    neighbor.days_infected = 0
                end
            end
        end

    elseif model.use_hospitalization && person.status == :H
        person.days_infected += 1

        # Check if agent recovers from hospital
        if person.days_infected >= model.days_to_hospital_recovery
            person.status = :R
        end
    end
end
