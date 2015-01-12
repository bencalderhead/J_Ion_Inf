using Distributions
using MCMCTypes

# Function for updating parameters - Gibbs step 1
include("MCMCUpdateParameters.jl")

# Function for updating indicator - Gibbs step 2
include("MCMCUpdateIndicator.jl")

# Function for the main MCMC routine
include("MCMCRun.jl")



# Define settings for simulation
OutputID            = "./Output/TestOutputFile2D"
Sampler             = "MH"
NumOfProposals      = 1
NumOfIterations     = 5000
InitialStepSize     = 1.0
ProposalCovariance  = [1.0 0.0; 0.0 1.0]
InitialiseFromPrior = false # Sample starting parameters from prior

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction     = UpdateParameters

# Define function for sampling the indicator variable (2nd Gibbs step)
SampleIndicatorFunction = SampleIndicator


############################
# Create an standard model #
############################

# Define the Model object
ModelType    = "Standard"
ModelName    = "Gaussian"
NumOfParas   = 2
ParaNames    = ["a" "b"]
DefaultParas = [0.0; 0.0] # Column vector

UsePrior     = false
Prior        = Array(PriorDistribution, NumOfParas)
Prior[1]     = PriorDistribution("Uniform", [0, 10])
Prior[2]     = PriorDistribution("Uniform", [0, 10])

# Specify LL function
function LL_Gaussian2d(x)
    # Input: must be a column 1d array
    # Output: Float64

    Mean = [0.0; 0.0] # Column vector
    Cov  = [1.0 0.0; 0.0 1.0]

    logpdf(MvNormal(Mean, Cov), x) # Returns a Float64
end

LLFun = LL_Gaussian2d


GaussModel = TargetOnlyModel( ModelType,
                              ModelName,
                              NumOfParas,
                              ParaNames,
                              DefaultParas,
                              UsePrior,
                              Prior,
                              LLFun)



############################
# Create simulation object #
############################

MySimulation = MCMCSimulation( OutputID,
                               Sampler,
                               NumOfProposals,
                               NumOfIterations,
                               InitialStepSize,
                               ProposalCovariance,
                               InitialiseFromPrior,
                               UpdateParasFunction,
                               SampleIndicatorFunction,
                               GaussModel )


# Run the MCMC code, passing in the MCMCSimulation object
MCMCRun( MySimulation )
