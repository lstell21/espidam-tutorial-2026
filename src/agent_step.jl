"""
    agent_step!(person::Person, model::AgentBasedModel)

Update the state of an agent in the ABM model for one time step.

# Arguments
- `person::Person`: The agent whose state needs to be updated.
- `model::AgentBasedModel`: The ABM model containing the agent.

# Description
- If `use_risk_switching` is enabled, susceptible HCWs switch to low-contact as the
  hospitalization rate rises (PPE adoption) and revert back to high-contact as it falls
  (HCWs let their guard down when hospitalizations decline).
- If infected, increment days infected and recover if enough time has passed.
- If infected, potentially transmit to nearby susceptibles; transmission probability is
  reduced by `low_risk_factor` when the infector has `contact_rate = :low`.
- If `use_hospitalization` is true, infected agents may be hospitalized on day 1.
  Hospitalized agents don't transmit and recover after `days_to_hospital_recovery` days.
"""
function agent_step!(person::Person, model::AgentBasedModel)
    # PPE adoption/reversal: high-contact susceptibles switch to low-contact
    # as hospitalization rate rises, and switch back as it falls (HCWs let guard down)
    if model.use_risk_switching && person.status == :S
        hosp_rate = model.hospitalized_count / model.n_nodes
        if person.contact_rate == :high
            # Rising hospitalizations → adopt PPE
            if rand(abmrng(model)) < hosp_rate * model.risk_switch_rate
                person.contact_rate = :low
            end
        else
            # Falling hospitalizations → revert to high-contact
            if rand(abmrng(model)) < (1 - hosp_rate) * model.risk_switch_rate
                person.contact_rate = :high
            end
        end
    end

    if person.status == :I
        person.days_infected += 1

        # Check if agent should be hospitalized (only on day 1 of infection)
        if model.use_hospitalization && person.days_infected == 1 && rand(abmrng(model)) < model.hospitalization_prob
            person.status = :H
            return
        end

        # Check if agent recovers from infection
        if person.days_infected >= model.days_to_recovered
            person.status = :R
            return
        end

        # Transmit to susceptible neighbors; transmission depends on infector's contact rate
        for neighbor in nearby_agents(person, model, 1)
            if neighbor.status == :S
                trans_prob = person.contact_rate == :high ? model.trans_prob : model.trans_prob * model.low_risk_factor
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
