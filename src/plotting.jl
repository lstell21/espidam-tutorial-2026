"""
    plot_degree_distribution(degree_distribution; network_type::Symbol)

Plot the degree distribution of a graph.

# Arguments
- `degree_distribution`: A dictionary mapping degrees to counts or DataFrame containing degree distribution
- `network_type`: Symbol indicating the type of network (for title)

# Returns
- `p`: A plot of the degree distribution.

# Example
```julia
plot_degree_distribution(graph_analysis["degree_distribution"], :random)
```
"""
function plot_degree_distribution(degree_distribution; network_type::Symbol)
    # Extract degrees and counts from different possible input formats
    if isa(degree_distribution, Dict)
        degrees = collect(keys(degree_distribution))
        counts = collect(values(degree_distribution))
    elseif isa(degree_distribution, DataFrame) && haskey(degree_distribution, :degree_distribution)
        dd = degree_distribution.degree_distribution[1]
        degrees = collect(keys(dd))
        counts = collect(values(dd))
    elseif isa(degree_distribution, DataFrame) && haskey(degree_distribution, :degree) && haskey(degree_distribution, :count)
        degrees = degree_distribution.degree
        counts = degree_distribution.count
    else
        error("Unsupported degree_distribution format")
    end
    
    # Create and return the plot
    return plot(degrees, counts, seriestype=:bar, 
                title="Degree Distribution ($(String(network_type)))", 
                xlabel="Degree", ylabel="Count", legend=:none)
end


"""
    plot_epidemic_trajectories(mdf, network_type::Symbol)

Plot the epidemic trajectories of susceptible, infected, and recovered individuals over time.

# Arguments
- `mdf`: A DataFrame containing the epidemic data with columns `:susceptible_count`, `:infected_count`, and `:recovered_count`.
- `network_type`: Symbol indicating the type of network (for title).

# Returns
- `p`: A plot of the epidemic trajectories.

# Example
```julia
plot_epidemic_trajectories(mdf, :random)
```
"""
function plot_epidemic_trajectories(mdf, network_type; title_suffix="")
    # Create the plot
    p = plot(mdf.time, mdf.susceptible_count, 
             label="Susceptible", 
             linewidth=2, 
             color=:blue,
             xlabel="Time (days)", 
             ylabel="Number of agents",
             title="SIHR Epidemic Dynamics - $(titlecase(string(network_type))) Network$(title_suffix)",
             legend=:right,
             size=(800, 500),
             margin=5mm)
    
    plot!(p, mdf.time, mdf.infected_count, 
          label="Infected", 
          linewidth=2, 
          color=:red)
    
    # Only plot hospitalized count if it exists in the dataframe
    if :hospitalized_count in names(mdf)
        plot!(p, mdf.time, mdf.hospitalized_count, 
              label="Hospitalized", 
              linewidth=2, 
              color=:orange)
    end
    
    plot!(p, mdf.time, mdf.recovered_count, 
          label="Recovered", 
          linewidth=2, 
          color=:green)
    
    return p
end


"""
    create_base_filename(model)

Create a consistent base filename from model parameters for saving files.
"""
function create_base_filename(model)
    return "$(model.network_type)_mdeg_$(model.mean_degree)_nn_$(model.n_nodes)_disp_$(model.dispersion)_pat0_$(model.patient_zero)_hicontact_$(model.high_contact)_hc_frac_$(model.fraction_high_contact)_lrf_$(model.low_risk_factor)_trans_$(model.trans_prob)"
end

"""
    plot_single_run(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000, 
                   dispersion::Float64=0.1, patient_zero::Symbol=:random, 
                   high_contact::Symbol=:random, fraction_high_contact::Float64=0.1, 
                   low_risk_factor::Float64=1.0, trans_prob::Float64=0.1, n_steps::Int=100, r̂=nothing, p̂=nothing)

Plot a single run of an epidemic simulation.

# Arguments
- `network_type::Symbol`: The type of network to use for the simulation. Possible values are `:random`, `:smallworld`, `:preferential`, `:configuration`, `:proportionatemixing`, `:edgelist`, or `:custom`.
- `mean_degree::Int`: The mean degree of the network. Default is 4.
- `n_nodes::Int`: The number of nodes in the network. Default is 1000.
- `dispersion::Float64`: The dispersion parameter for the network. Default is 0.1.
- `patient_zero::Symbol`: The type of patient zero to use for the simulation. Default is `:random`.
- `high_contact::Symbol`: How high-risk individuals are distributed. Default is `:random`.
- `fraction_high_contact::Float64`: The fraction of high-risk individuals in the population. Default is 1.0.
- `low_risk_factor::Float64`: The factor by which low-risk individuals' transmission probability is multiplied. Must be between 0 and 1. Default is 1.0.
- `trans_prob::Float64`: The transmission probability. Default is 0.1.
- `n_steps::Int`: The number of simulation steps to run. Default is 100.
- `r̂`: The r parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`. Default is nothing.
- `p̂`: The p parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`. Default is nothing.

# Returns
- `plotdynamics`: A plot of the epidemic trajectories.
- `plotdegdist`: A plot of the degree distribution.
- `combined_plot`: A combined plot with both plots side-by-side.

# Example
```julia
dynamics_plot, degdist_plot, combined_plot = plot_single_run(network_type=:random, mean_degree=2)
```
"""
function plot_single_run(; network_type::Symbol, mean_degree::Int=4, n_nodes::Int=1000,
                       dispersion::Float64=0.1, patient_zero::Symbol=:random,
                       high_contact::Symbol=:random, fraction_high_contact::Float64=1.0,
                       low_risk_factor::Float64=1.0, trans_prob::Float64=0.1, n_steps::Int=100,
                       r̂=nothing, p̂=nothing, use_hospitalization::Bool=true,
                       figures_dir::String="figures")
    # Validate low_risk_factor
    if !(0 <= low_risk_factor <= 1)
        error("low_risk_factor must be between 0 and 1, got $low_risk_factor")
    end

    # Initialize model and run simulation
    model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero,
                     high_contact, fraction_high_contact, low_risk_factor, trans_prob, r̂, p̂,
                     use_hospitalization)

    # Define adata and mdata locally to avoid relying on global variables
    adata = [:status]
    mdata = use_hospitalization ?
        [:susceptible_count, :infected_count, :hospitalized_count, :recovered_count] :
        [:susceptible_count, :infected_count, :recovered_count]
    
    # Run the simulation
    _, mdf = run!(model, n_steps; adata, mdata)
    
    # Get graph analysis results
    graph_analysis = analyze_graph(model.graph)
    
    # Create base filename for consistency
    base_filename = create_base_filename(model)
    
    # Create plots
    plotdynamics = plot_epidemic_trajectories(mdf, model.network_type)
    plotdegdist = plot_degree_distribution(graph_analysis["degree_distribution"]; network_type=model.network_type)
    
    # Create a combined plot for side-by-side visualization
    combined_plot = plot(plotdynamics, plotdegdist, layout=(1,2), size=(1000, 400), margin=5mm,
                       title=["Epidemic Dynamics ($(String(model.network_type)))" "Degree Distribution ($(String(model.network_type)))"])
    mkpath(figures_dir)
    savefig(combined_plot, joinpath(figures_dir, "combined_plot_$(base_filename).pdf"))
    
    # Return all three plots
    return plotdynamics, plotdegdist, combined_plot
end


"""
    create_comparison_box_plot(formatted_data, title, ylabel; filename=nothing)

Create a standardized box plot for comparison data.

# Arguments
- `formatted_data`: A DataFrame with columns `value` and `network_type`
- `title`: Plot title
- `ylabel`: Y-axis label
- `filename`: Optional filename to save the plot

# Returns
- A boxplot comparing the data across network types
"""
function create_comparison_box_plot(formatted_data, title, ylabel; filename=nothing, colors=nothing, network_order=nothing)
    # Convert network_type to a categorical variable with ordered levels if order is provided
    if !isnothing(network_order)
        formatted_data.network_type = CategoricalArray(formatted_data.network_type, 
                                                      ordered=true, 
                                                      levels=network_order)
    end
    
    # If colors are provided, use them for the boxplot
    p = boxplot(
        formatted_data.network_type, formatted_data.value,
        title=title,
        xlabel="Network Type",
        ylabel=ylabel,
        legend=false,
        size=(600, 300),
        margin=5mm,
        guidefontsize=9,
        titlefontsize=10,
        xtickfontsize=8,
        xrotation=30,
        linewidth=1,
        left_margin=8mm,
        outliers=false
        )
    
    if !isnothing(filename)
        savefig(p, filename)
    end
    
    return p
end


"""
    run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4, 
                          n_nodes::Int=1000, dispersion::Float64=0.1, 
                          patient_zero::Symbol=:random, high_contact::Symbol=:random, 
                          fraction_high_contact::Float64=1.0, low_risk_factor::Float64=1.0,
                          trans_prob::Float64=0.1, n_steps::Int=100, boxplot_colors=nothing, r̂=nothing, p̂=nothing)

Run simulations for multiple network types and generate comparison plots.
"""
function run_and_plot_comparison(; network_types::Vector{Symbol}, mean_degree::Int=4,
                               n_nodes::Int=1000, dispersion::Float64=0.1,
                               patient_zero::Symbol=:random, high_contact::Symbol=:random,
                               fraction_high_contact::Float64=1.0, low_risk_factor::Float64=1.0,
                               trans_prob::Float64=0.1, n_steps::Int=100, boxplot_colors=nothing,
                               r̂=nothing, p̂=nothing, use_hospitalization::Bool=true,
                               figures_dir::String="figures", data_dir::String="data",
                               output_dir_path::String="output", save_data::Bool=false)
    # Validate low_risk_factor
    if !(0 <= low_risk_factor <= 1)
        error("low_risk_factor must be between 0 and 1, got $low_risk_factor")
    end
    
    # Store results for each network type
    all_results = Dict()
    
    # Prepare empty DataFrames for each metric
    duration_data = DataFrame(value=Float64[], network_type=String[])
    max_infected_data = DataFrame(value=Float64[], network_type=String[])
    susceptible_remaining_data = DataFrame(value=Float64[], network_type=String[])
    
    # Run simulations for each network type
    for network_type in network_types
        println("Running simulations for $(network_type) network...")
        model = initialize(; network_type, mean_degree, n_nodes, dispersion, patient_zero, 
                          high_contact, fraction_high_contact, low_risk_factor, trans_prob, r̂, p̂)
        
        # Run simulations
        multiple_runs = run_simulations(; network_type, mean_degree, n_nodes, dispersion, 
                                       patient_zero, high_contact, fraction_high_contact, 
                                       low_risk_factor, trans_prob, n_steps, r̂, p̂, use_hospitalization)
        
        # Process results
        grouped_data = groupby(multiple_runs, [:seed])
        
        # Calculate duration metrics properly
        final_results = combine(grouped_data) do df
            # Find the first step where infection appears (should be step 1 normally)
            first_infected_step = findfirst(df.infected_count .> 0)
            
            # Find the last step where infection is present
            last_infected_step = findlast(df.infected_count .> 0)
            
            # Calculate the duration - if epidemic never ends, use total steps
            duration = if isnothing(last_infected_step)
                length(df.infected_count)  # Epidemic didn't end
            else
                last_infected_step - first_infected_step + 1  # +1 to include both endpoints
            end
            
            # Get maximum infected count
            max_infected = maximum(df.infected_count)
            
            return DataFrame(
                first_to_last_infected = duration,
                max_infected = max_infected
            )
        end
        
        last_rows = combine(grouped_data, names(multiple_runs) .=> last)
        final_results[!, :susceptible_fraction_remaining] = last_rows.susceptible_count_last ./ 
            (last_rows.susceptible_count_last + last_rows.infected_count_last + last_rows.recovered_count_last)
        
        if save_data
            base_filename = "$(network_type)_mdeg_$(mean_degree)_nn_$(n_nodes)_disp_$(dispersion)_pat0_$(patient_zero)_hicontact_$(high_contact)_hc_frac_$(fraction_high_contact)_lrf_$(low_risk_factor)_trans_$(trans_prob)"
            mkpath(data_dir)
            mkpath(output_dir_path)
            CSV.write(joinpath(data_dir, "simulation_results_$(base_filename).csv"), multiple_runs)
            CSV.write(joinpath(output_dir_path, "final_results_$(base_filename).csv"), final_results)
        end
        
        # Store results for box plots
        all_results[network_type] = final_results
        
        # Add data for each metric to the appropriate DataFrame with network type as a label
        network_label = titlecase(String(network_type))
        for value in final_results.first_to_last_infected
            push!(duration_data, (value, network_label))
        end
        
        for value in final_results.max_infected
            push!(max_infected_data, (value, network_label))
        end
        
        for value in final_results.susceptible_fraction_remaining
            push!(susceptible_remaining_data, (value, network_label))
        end
    end
    
    # Create an ordered list of network types for consistent ordering across plots
    network_order = [titlecase(String(nt)) for nt in network_types]
    
    # Generate comparison plots using the pre-formatted data and colors
    mkpath(figures_dir)
    net_types_str = join([String(nt) for nt in network_types], "_")

    duration_comparison = create_comparison_box_plot(
        duration_data,
        "Epidemic Duration Comparison\n(mean degree: $(mean_degree))",
        "Duration (steps)",
        filename=joinpath(figures_dir, "duration_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf"),
        colors=boxplot_colors,
        network_order=network_order
    )

    max_infected_comparison = create_comparison_box_plot(
        max_infected_data,
        "Maximum Infected Comparison\n(mean degree: $(mean_degree))",
        "Maximum Number of Infected",
        filename=joinpath(figures_dir, "max_infected_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf"),
        colors=boxplot_colors,
        network_order=network_order
    )

    sfr_comparison = create_comparison_box_plot(
        susceptible_remaining_data,
        "Susceptible Fraction Remaining Comparison\n(mean degree: $(mean_degree))",
        "Susceptible Fraction Remaining",
        filename=joinpath(figures_dir, "sfr_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf"),
        colors=boxplot_colors,
        network_order=network_order
    )

    combined_comparison = plot(
        duration_comparison, max_infected_comparison, sfr_comparison,
        layout=(1,3),
        size=(1200, 500),
        margin=3mm,
        bottom_margin=10mm,
        top_margin=5mm,
        title=["Epidemic Duration" "Maximum Infected" "Susceptible Fraction Remaining"],
        titlefontsize=10,
        plot_title="Epidemic Comparison (mean degree: $(mean_degree))",
        plot_titlefontsize=12,
        left_margin=8mm,
    )
    savefig(combined_comparison, joinpath(figures_dir, "combined_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf"))
    
    return combined_comparison
end


"""
    plot_centrality_comparison(;network_types=[:random, :smallworld, :preferential], 
                              mean_degree=4, n_nodes=1000, link_axes=false, r̂=nothing, p̂=nothing)

Plot boxplots comparing centrality measures across different network types.

# Arguments
- `network_types`: Vector of symbols representing the network types to compare (any number supported)
- `mean_degree`: Mean degree for network generation
- `n_nodes`: Number of nodes in each network
- `link_axes`: Boolean indicating whether to link y-axes across plots for easier comparison (default: false)
- `r̂`: The r parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`
- `p̂`: The p parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`

# Returns
- A combined plot showing boxplots of centrality measures for each network type

# Example
```julia
# Default visualization with independent y-axes for three network types
centrality_comparison = plot_centrality_comparison()

# With linked y-axes for direct comparison with two network types
centrality_comparison = plot_centrality_comparison(
    network_types=[:random, :preferential],
    link_axes=true
)

# With five different network types
centrality_comparison = plot_centrality_comparison(
    network_types=[:random, :smallworld, :preferential, :configuration, :proportionatemixing]
)
```
"""
function plot_centrality_comparison(;network_types=[:random, :smallworld, :preferential],
                                   mean_degree=4, n_nodes=1000, link_axes=false,
                                   r̂=nothing, p̂=nothing, figures_dir::String="figures")
    # Initialize empty DataFrames to store the centrality data
    centrality_data = Dict()
    
    # Generate and analyze each network type
    for nt in network_types
        model = initialize(; network_type=nt, mean_degree=mean_degree, n_nodes=n_nodes, r̂=r̂, p̂=p̂)
        analysis = analyze_graph(model.graph)
        centrality_data[nt] = analysis["centrality"]
    end
    
    # Determine the optimal plot layout based on the number of network types
    n_types = length(network_types)
    
    if n_types == 1
        # For a single network type, use 1x1
        plot_layout = (1, 1)
        plot_width = 600
    elseif n_types == 2
        # For two network types, use 1x2
        plot_layout = (1, 2)
        plot_width = 1000
    elseif n_types <= 4
        # For 3-4 network types, use 1xN
        plot_layout = (1, n_types)
        plot_width = min(1600, 500 * n_types)
    else
        # For more than 4 network types, use a more compact grid layout
        n_cols = ceil(Int, sqrt(n_types))
        n_rows = ceil(Int, n_types / n_cols)
        plot_layout = (n_rows, n_cols)
        plot_width = min(1800, 450 * n_cols)
    end
    
    # Fixed plot height per row
    plot_height_per_row = 450
    plot_height = plot_height_per_row * first(plot_layout)
    
    # Create plots for each centrality measure
    plots = []
    
    # Calculate y-axis limits for each network type independently
    # This allows better visualization of the specific distributions
    y_limits = Dict()
    
    for nt in network_types
        df = centrality_data[nt]
        
        # Calculate quartiles for each centrality measure
        q1_degree = quantile(df.degree_centrality, 0.25)
        q3_degree = quantile(df.degree_centrality, 0.75)
        iqr_degree = q3_degree - q1_degree
        upper_whisker_degree = min(maximum(df.degree_centrality), q3_degree + 1.5 * iqr_degree)
        
        q1_betweenness = quantile(df.betweenness_centrality, 0.25)
        q3_betweenness = quantile(df.betweenness_centrality, 0.75)
        iqr_betweenness = q3_betweenness - q1_betweenness
        upper_whisker_betweenness = min(maximum(df.betweenness_centrality), q3_betweenness + 1.5 * iqr_betweenness)
        
        q1_closeness = quantile(df.closeness_centrality, 0.25)
        q3_closeness = quantile(df.closeness_centrality, 0.75)
        iqr_closeness = q3_closeness - q1_closeness
        upper_whisker_closeness = min(maximum(df.closeness_centrality), q3_closeness + 1.5 * iqr_closeness)
        
        q1_eigenvector = quantile(df.eigenvector_centrality, 0.25)
        q3_eigenvector = quantile(df.eigenvector_centrality, 0.75)
        iqr_eigenvector = q3_eigenvector - q1_eigenvector
        upper_whisker_eigenvector = min(maximum(df.eigenvector_centrality), q3_eigenvector + 1.5 * iqr_eigenvector)
        
        # Use upper whiskers (excluding outliers) with a small buffer for y-axis limits
        degree_max = upper_whisker_degree * 1.15
        betweenness_max = upper_whisker_betweenness * 1.15
        closeness_max = upper_whisker_closeness * 1.15
        eigenvector_max = upper_whisker_eigenvector * 1.15
        
        # Store the maximum overall for this network type
        y_limits[nt] = max(degree_max, betweenness_max, closeness_max, eigenvector_max)
    end
    
    # If link_axes is true, find the global maximum across all networks
    global_y_max = link_axes ? maximum(values(y_limits)) : 0
    
    # Create a DataFrame for each network type to use with StatsPlots groupedboxplot
    for (i, nt) in enumerate(network_types)
        df = centrality_data[nt]
        network_name = titlecase(String(nt))
        
        # Create a DataFrame for boxplot data
        boxplot_data = DataFrame(
            centrality_type = vcat(
                fill("Degree", nrow(df)),
                fill("Betweenness", nrow(df)),
                fill("Closeness", nrow(df)),
                fill("Eigenvector", nrow(df))
            ),
            value = vcat(
                df.degree_centrality,
                df.betweenness_centrality,
                df.closeness_centrality,
                df.eigenvector_centrality
            )
        )
        
        # Convert centrality_type to a categorical array with ordered levels
        boxplot_data.centrality_type = CategoricalArray(
            boxplot_data.centrality_type, 
            ordered=true, 
            levels=["Degree", "Betweenness", "Closeness", "Eigenvector"]
        )
        
        # Set y-axis limits based on link_axes parameter
        y_max = link_axes ? global_y_max : y_limits[nt]
        
        # Create a more precisely aligned boxplot using @df macro
        p = @df boxplot_data boxplot(
            :centrality_type, 
            :value,
            title=network_name,
            legend=false,
            outliers=false,  # Hide outliers to improve readability
            marker=(0.5, :circle, 0.3),
            alpha=0.7,       
            guidefontsize=9,
            titlefontsize=10,
            widen=true,      
            bar_width=0.7,   
            xticks=([1,2,3,4], ["Degree", "Betweenness", "Closeness", "Eigenvector"]),
            xrotation=30,    
            ylims=(0, y_max),
            ylabel="Centrality Value",
            dpi=300,
            margin=8mm,      # Add individual plot margins for better spacing
            bottom_margin=10mm,  # Add extra space at the bottom for x-axis labels
            left_margin=8mm      # Add extra space for y-axis labels
        )
        
        push!(plots, p)
    end
    
    # Handle case where we might have more layout slots than plots
    if length(plots) < prod(plot_layout)
        # Fill remaining slots with empty plots
        for i in length(plots)+1:prod(plot_layout)
            push!(plots, plot(framestyle=:none, grid=false, showaxis=false))
        end
    end
    
    # Combine plots in a single figure with improved margins
    combined_plot = plot(plots..., 
                      layout=plot_layout, 
                      size=(plot_width, plot_height),
                      bottom_margin=12mm, # Increase bottom margin for x-axis labels
                      top_margin=10mm,    # Add more top margin for titles
                      xtickfontsize=9,
                      link=link_axes ? :y : :none,
                      title_position=:center)
    
    net_types_str = join([String(nt) for nt in network_types], "_")
    mkpath(figures_dir)
    savefig(combined_plot, joinpath(figures_dir, "centrality_comparison_$(net_types_str)_mdeg_$(mean_degree).pdf"))
    
    return combined_plot
end


"""
    plot_network_metrics_comparison(;network_types=[:random, :smallworld, :preferential], 
                                  mean_degree=4, n_nodes=1000, r̂=nothing, p̂=nothing)

Create a horizontal bar plot comparing structural metrics across different network types.

# Arguments
- `network_types`: Vector of symbols representing the network types to compare
- `mean_degree`: Mean degree for network generation
- `n_nodes`: Number of nodes in each network
- `r̂`: The r parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`
- `p̂`: The p parameter for negative binomial distribution, used only when `network_type` is `:proportionatemixing`

# Returns
- A horizontal bar chart comparing network metrics across network types

# Example
```julia
# Default comparison of three network types
metrics_plot = plot_network_metrics_comparison()

# Custom comparison with two network types
metrics_plot = plot_network_metrics_comparison(
    network_types=[:random, :preferential],
    mean_degree=6
)
```
"""
function plot_network_metrics_comparison(;network_types=[:random, :smallworld, :preferential],
                                       mean_degree=4, n_nodes=1000, r̂=nothing, p̂=nothing,
                                       figures_dir::String="figures")
    # Define a colorblind-friendly palette for the five network types
    # Using a modified version of Wong's palette which is colorblind-friendly
    network_color_map = Dict(
        :random => RGB(0/255, 114/255, 178/255),        # Blue
        :smallworld => RGB(230/255, 159/255, 0/255),    # Orange
        :preferential => RGB(0/255, 158/255, 115/255),  # Green
        :configuration => RGB(204/255, 121/255, 167/255), # Purple
        :proportionatemixing => RGB(213/255, 94/255, 0/255), # Red-orange
        :proportionate => RGB(213/255, 94/255, 0/255),   # Same as proportionatemixing (alternative name)
        :edgelist => RGB(86/255, 180/255, 233/255),     # Light blue
        :custom => RGB(240/255, 228/255, 66/255)        # Yellow
    )
    
    # Create storage for metrics
    density_values = Float64[]
    clustering_values = Float64[]
    assortativity_values = Float64[]
    degree_cent_values = Float64[]
    betweenness_cent_values = Float64[]
    closeness_cent_values = Float64[]
    eigenvector_cent_values = Float64[]
    network_labels = String[]
    network_colors = []

    # Collect metrics for each network type
    for nt in network_types
        println("Analyzing $(nt) network...")
        model = initialize(; network_type=nt, mean_degree=mean_degree, n_nodes=n_nodes, r̂=r̂, p̂=p̂)
        analysis = analyze_graph(model.graph)
        
        # Extract key metrics
        density = analysis["summary"][analysis["summary"].metric .== "Density", :value][1]
        clustering = analysis["summary"][analysis["summary"].metric .== "Clustering Coefficient", :value][1]
        assortativity = analysis["summary"][analysis["summary"].metric .== "Assortativity", :value][1]
        degree_cent = mean(analysis["centrality"].degree_centrality)
        betweenness_cent = mean(analysis["centrality"].betweenness_centrality)
        closeness_cent = mean(analysis["centrality"].closeness_centrality)
        eigenvector_cent = mean(analysis["centrality"].eigenvector_centrality)
        
        # Store values
        push!(density_values, density)
        push!(clustering_values, clustering)
        push!(assortativity_values, assortativity)
        push!(degree_cent_values, degree_cent)
        push!(betweenness_cent_values, betweenness_cent)
        push!(closeness_cent_values, closeness_cent)
        push!(eigenvector_cent_values, eigenvector_cent)
        
        # Get consistent color for this network type
        push!(network_colors, get(network_color_map, nt, :gray))
        
        # Create display label
        label = if nt == :random
            "Random"
        elseif nt == :smallworld
            "Small-World"
        elseif nt == :preferential
            "Preferential Attachment"
        elseif nt == :configuration
            "Configuration"
        elseif nt == :proportionatemixing || nt == :proportionate
            "Proportionate Mixing"
        elseif nt == :edgelist
            "Edgelist"
        elseif nt == :custom
            "Custom"
        else
            String(nt)
        end
        push!(network_labels, label)
    end

    # Create DataFrame with metrics
    df = DataFrame(
        "metric" => [
            "Density", 
            "Clustering Coefficient", 
            "Assortativity", 
            "Degree Centrality", 
            "Betweenness Centrality", 
            "Closeness Centrality", 
            "Eigenvector Centrality"
        ],
        "position" => 1:7
    )

    # Add columns for each network type, maintaining the order from the network_types argument
    for (i, label) in enumerate(network_labels)
        df[!, label] = [
            density_values[i],
            clustering_values[i],
            assortativity_values[i],
            degree_cent_values[i],
            betweenness_cent_values[i],
            closeness_cent_values[i],
            eigenvector_cent_values[i]
        ]
    end

    # Reshape for plotting
    df_long = stack(df, Not(["metric", "position"]), variable_name = "Model", value_name = "Value")
    
    # Convert Model to a categorical variable with ordered levels to preserve the order
    # from the network_types argument
    df_long.Model = CategoricalArray(df_long.Model, ordered=true, levels=network_labels)

    # Create grouped bar plot with improved spacing and margins
    metrics_plot = @df df_long groupedbar(
        :metric, 
        :Value,
        group = :Model,
        bar_position = :dodge,
        xlabel = "Metric",
        ylabel = "Value",
        title = "Network Metrics by Model Type",
        legend = :topright,
        size = (900, 550),
        ylims = (-0.15, 0.6),
        left_margin = 5mm,
        top_margin = 5mm,
        bottom_margin = 15mm,
        xrotation = 45,
        palette = network_colors  # Use our defined colors
    )

    network_types_str = join([String(nt) for nt in network_types], "_")
    mkpath(figures_dir)
    savefig(metrics_plot, joinpath(figures_dir, "network_metrics_comparison_$(network_types_str)_mdeg_$(mean_degree).pdf"))

    return metrics_plot
end

"""
    plot_epidemic_comparison(results_dict; save_plots=true)

Plot the comparison of epidemic trajectories across different network types.

# Arguments
- `results_dict`: A dictionary containing the results of the epidemic simulations for different network types.
- `save_plots`: Boolean indicating whether to save the plot as a PDF file. Default is true.

# Returns
- `p`: A plot object representing the epidemic comparison plot.

# Example
```julia
plot_epidemic_comparison(results_dict)
```
"""
function plot_epidemic_comparison(results_dict; save_plots=true)
    # Create comparison plot
    p = plot(xlabel="Time (days)", 
             ylabel="Number of agents",
             title="SIHR Epidemic Comparison Across Network Types",
             legend=:outertopright,
             size=(1000, 600),
             margin=8mm)
    
    colors = Dict(:random => :blue, :smallworld => :red, :preferential => :green, 
                  :proportionate => :purple, :edgelist => :orange)
    
    for (network_type, result) in results_dict
        if haskey(result, "mean_df")
            mdf = result["mean_df"]
            
            # Plot mean trajectories with confidence intervals if available
            plot!(p, mdf.time, mdf.infected_count_mean, 
                  label="$(titlecase(string(network_type))) - Infected",
                  linewidth=2, 
                  color=colors[network_type],
                  linestyle=:solid)
            
            plot!(p, mdf.time, mdf.hospitalized_count_mean, 
                  label="$(titlecase(string(network_type))) - Hospitalized",
                  linewidth=2, 
                  color=colors[network_type],
                  linestyle=:dash)
            
            # Add confidence intervals if available
            if haskey(mdf, :infected_count_ci_lower)
                plot!(p, mdf.time, mdf.infected_count_ci_lower, 
                      fillto=mdf.infected_count_ci_upper,
                      alpha=0.2, 
                      color=colors[network_type],
                      label="")
            end
        end
    end
    
    # Save plot if requested
    if save_plots
        savefig(p, "figures/epidemic_comparison_SIHR.pdf")
        println("Epidemic comparison plot saved to figures/epidemic_comparison_SIHR.pdf")
    end
    
    return p
end