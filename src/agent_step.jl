"""
    agent_step!(person::Person, model::AgentBasedModel)

Update the state of an agent in the ABM model for one time step.

# Arguments
- `person::Person`: The agent whose state needs to be updated.
- `model::AgentBasedModel`: The ABM model containing the agent.

# Description
- If infected, increment days infected and recover if enough time has passed.
- If infected, potentially transmit to nearby susceptibles. The base per-contact
  transmission probability is `trans_prob` for high-contact infectors (HCWs) and
  `trans_prob * low_risk_factor` for low-contact infectors (patients).
- If `use_ppe_adoption` is enabled, HCW transmission is additionally scaled down by a
  behavioral "awareness" factor that grows with the observed hospitalization burden
  (PPE adoption). The factor is a saturating function of the hospitalized fraction, so
  it rises as hospitalizations climb and relaxes back as they fall — no per-agent state
  is switched, only the effective transmission probability is modulated.
- If `use_hospitalization` is true, infected agents may be hospitalized on day 1.
  Hospitalized agents don't transmit and recover after `days_to_hospital_recovery` days.
"""
function agent_step!(person::Person, model::AgentBasedModel)
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

        # Effective per-contact transmission probability for this infector.
        # Identity gap (HCW vs patient) is fixed via low_risk_factor; PPE adoption is a
        # dynamic, severity-driven reduction applied to HCW transmission only.
        trans_prob = person.contact_rate == :high ? model.trans_prob : model.trans_prob * model.low_risk_factor
        if model.use_ppe_adoption && person.contact_rate == :high
            signal = model.hospitalized_count / model.n_nodes
            reduction = model.max_ppe_reduction * signal / (signal + model.ppe_half_saturation)
            trans_prob *= (1.0 - reduction)
        end

        # Transmit to susceptible neighbors
        for neighbor in nearby_agents(person, model, 1)
            if neighbor.status == :S
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