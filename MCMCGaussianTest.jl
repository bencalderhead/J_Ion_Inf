using Distributions
using MCMCTypes

# Function for updating parameters - Gibbs step 1
include("MCMCUpdateParameters.jl")

# Function for updating indicator - Gibbs step 2
include("MCMCUpdateIndicator.jl")

# Function for the main MCMC routine
include("MCMCRun.jl")



# Define settings for simulation
OutputID            = "./Output/TestOutputFile_1D"
Sampler             = "MH"
NumOfProposals      = 1
NumOfIterations     = 5000
InitialStepSize     = 1.0
ProposalCovariance  = [1.0]
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
NumOfParas   = 1
ParaNames    = ["a"]
DefaultParas = [0.0]

UsePrior     = false
Prior        = Array(PriorDistribution, 1)
Prior[1]     = PriorDistribution("Uniform", [0, 10])


# Specify LL function
function LL_Gaussian(x, Model)
    Mean = 0.0
    Std  = 1.0

    logpdf(Normal(Mean, Std), x[1]) # Note that this univariate distribution returns a Float if x is a Float
end


LLFun = LL_Gaussian



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
