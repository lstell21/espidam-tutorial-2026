"""
    agent_step!(person::Person, model::AgentBasedModel)

Update the state of an agent in the ABM model for one time step.

# Arguments
- `person::Person`: The agent whose state needs to be updated.
- `model::AgentBasedModel`: The ABM model containing the agent.

# Description
- If `use_behavior_adoption` is enabled, every agent that has not yet adopted
  contact-restricting behavior draws (once per step) whether to adopt it. The adoption
  probability is a logistic function of the observed hospitalized fraction
  `s = hospitalized_count / n_nodes`:
  `behavior_max_prob / (1 + exp(-behavior_alpha * (s - behavior_beta)))`.
  Adoption is absorbing (once adopted, it stays), creating a negative feedback loop:
  more hospitalizations → more adoption → fewer contacts → fewer infections.
- If infected, increment days infected and recover if enough time has passed.
- If infected, potentially transmit to nearby susceptibles. The base per-contact
  transmission probability is `trans_prob` for high-contact infectors (HCWs) and
  `trans_prob * low_risk_factor` for low-contact infectors (patients).
- A contact involves two agents, so each endpoint that restricts its contacts scales the
  transmission probability for that contact by `(1 - contact_reduction)`. This protects
  susceptible adopters as well as throttling infectious ones.
- If `use_hospitalization` is true, infected agents may be hospitalized on day
  `model.hospitalization_day` of infection (default 1). Hospitalized agents
  don't transmit and recover after `days_to_hospital_recovery` days.
"""
function agent_step!(person::Person, model::AgentBasedModel)
    # Hospitalization-driven adoption of contact-restricting behavior (absorbing).
    # Applies to every agent regardless of status; the hospitalized fraction is the
    # value recorded at the end of the previous step, so all agents see the same signal.
    if model.use_behavior_adoption && !model.restricts_contact[person.id]
        s = model.hospitalized_count / model.n_nodes
        adopt_prob = model.behavior_max_prob / (1.0 + exp(-model.behavior_alpha * (s - model.behavior_beta)))
        if rand(abmrng(model)) < adopt_prob
            model.restricts_contact[person.id] = true
        end
    end

    if person.status == :I
        person.days_infected += 1

        # Check if agent should be hospitalized (only on the configured day of infection)
        if model.use_hospitalization && person.days_infected == model.hospitalization_day && rand(abmrng(model)) < model.hospitalization_prob
            person.status = :H
            return
        end

        # Check if agent recovers from infection
        if person.days_infected >= model.days_to_recovered
            person.status = :R
            return
        end

        # Effective per-contact transmission probability for this infector. The identity
        # gap (HCW vs patient) is fixed via low_risk_factor; an infector that restricts
        # its contacts scales this down by (1 - contact_reduction).
        trans_prob = person.contact_rate == :high ? model.trans_prob : model.trans_prob * model.low_risk_factor
        if model.use_behavior_adoption && model.restricts_contact[person.id]
            trans_prob *= (1.0 - model.contact_reduction)
        end

        # Transmit to susceptible neighbors. A susceptible neighbor that restricts its
        # own contacts further reduces the chance of transmission across that contact.
        for neighbor in nearby_agents(person, model, 1)
            if neighbor.status == :S
                edge_prob = trans_prob
                if model.use_behavior_adoption && model.restricts_contact[neighbor.id]
                    edge_prob *= (1.0 - model.contact_reduction)
                end
                if rand(abmrng(model)) < edge_prob
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