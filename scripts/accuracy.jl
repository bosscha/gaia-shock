## Simulations to measure accuracy of the method with parallax noise
##
#### Gemini produced ...
using Distributions , Random
using Statistics

function simulate_parallax_noise(true_parallax::Float64, gaia_mag::Float64;
                                 num_observations::Int=100,
                                 base_parallax_error::Float64=0.03,  # Base error in milliarcseconds (mas)
                                 mag_dependence_factor::Float64=0.1, # Controls how much magnitude impacts error
                                 scan_angle_effect::Bool=true,        # Include scan-angle-dependent error
                                 scan_angle_std::Float64=0.015,    # Standard deviation of scan angle effect (mas)
                                 systematic_error_floor::Float64 = 0.01, # Minimum error (mas), even for bright stars
                                 add_outliers::Bool=false,          # Simulate occasional large errors
                                 outlier_probability::Float64=0.01,  # Probability of an outlier
                                 outlier_std_factor::Float64=10.0,   # Outlier standard deviation multiplier
                                 rng::AbstractRNG = Random.GLOBAL_RNG #allow passing a random number generator
                                 )

    """
    Simulates parallax measurements with noise, mimicking Gaia data characteristics.

    Args:
        true_parallax: The true parallax of the star in milliarcseconds (mas).
        gaia_mag: The Gaia G-band magnitude of the star.
        num_observations: The number of simulated parallax observations.
        base_parallax_error: The base parallax error for a bright star (mas).
        mag_dependence_factor:  A factor controlling the magnitude dependence of the error.
                                Higher values mean magnitude has a larger impact.
        scan_angle_effect: Whether to include a scan-angle-dependent error component.
        scan_angle_std: Standard deviation of the scan angle effect (mas).
        systematic_error_floor: Minimum error achievable, even for very bright stars (mas).
        add_outliers:  Whether to include occasional large outlier measurements.
        outlier_probability: Probability of a single observation being an outlier.
        outlier_std_factor:  How many times larger the outlier standard deviation is
                              compared to the normal error.
        rng:  Random number generator.  Defaults to the global RNG.  Pass a specific
             instance (e.g., `MersenneTwister(1234)`) for reproducible results.


    Returns:
        A vector of simulated parallax measurements (mas).  The length of the vector
        will be equal to `num_observations`.

    Example:
    ```julia
    simulated_parallaxes = simulate_parallax_noise(1.0, 15.0)  # True parallax = 1 mas, G=15 mag
    println("Mean simulated parallax: ", mean(simulated_parallaxes))
    println("Standard deviation: ", std(simulated_parallaxes))

    # Reproducible example with a specific seed:
    using Random
    rng = MersenneTwister(42);  # Use a specific seed
    simulated_parallaxes_reproducible = simulate_parallax_noise(1.0, 15.0; rng=rng)
    ```
    """

    # Magnitude-dependent error.  Brighter stars have smaller errors.
    mag_error = base_parallax_error * (1.0 + mag_dependence_factor * (gaia_mag - 12.0))
    mag_error = max(mag_error, systematic_error_floor) # Apply error floor


    # Generate the simulated measurements
    measurements = Vector{Float64}(undef, num_observations)
    for i in 1:num_observations
        # Basic parallax error (normally distributed)
        parallax_error = rand(rng, Normal(0.0, mag_error))

        # Scan angle effect (if enabled)
        if scan_angle_effect
            scan_angle_error = rand(rng, Normal(0.0, scan_angle_std))
            parallax_error += scan_angle_error
        end

        # Add outliers
        if add_outliers && rand(rng) < outlier_probability
             parallax_error += rand(rng, Normal(0.0, mag_error * outlier_std_factor))
        end

        measurements[i] = true_parallax + parallax_error
    end

    return measurements
end

### end Gemini produced...




#################################### MAIN###############################
let
    args = parse_commandline()
    debug_red("# Gemini produced of parallax errors ....")