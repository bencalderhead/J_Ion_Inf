using Distributions
using MCMCTypes

# Function for updating parameters - Gibbs step 1
include("MCMCUpdateParameters.jl")

# Function for updating indicator - Gibbs step 2
include("MCMCUpdateIndicator.jl")

# Function for the main MCMC routine
include("MCMCRun.jl")



# Define settings for simulation
OutputID            = "./Output/BivariateGibbsOutput"
Sampler             = "Gibbs"
NumOfProposals      = 1
NumOfIterations     = 1000000
InitialStepSize     = 1.0
ProposalCovariance  = []
InitialiseFromPrior = false # Sample starting parameters from prior

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction     = UpdateParameters

# Define function for sampling the indicator variable (2nd Gibbs step)
SampleIndicatorFunction = SampleIndicator


############################
# Create an standard model #
############################

# Define the Model object
ModelType    = "Bivariate"
ModelName    = "Gaussian"
NumOfParas   = 2
ParaNames    = ["a" "b"]
DefaultParas = [0.0; 0.0] # Column vector

UsePrior     = false
Prior        = Array(Distribution, NumOfParas)
Prior[1]     = Uniform(0, 10)
Prior[2]     = Uniform(0, 10)

# Specify LL function
function LL_Gaussian2d(x, Model::GaussianBivariate)
    # Input: must be a column 1d array
    # Output: Float64

    logpdf(MvNormal([Model.TargetMean1; Model.TargetMean2], Model.TargetCov), x) # Returns a Float64
end

LLFun = LL_Gaussian2d

TargetMean1 = 0.0
TargetMean2 = 0.0
TargetStd1  = 1.0
TargetStd2  = 1.0
Rho         = 0.99


GaussModel = GaussianBivariate( ModelType,
                                ModelName,
                                NumOfParas,
                                ParaNames,
                                DefaultParas,
                                UsePrior,
                                Prior,
                                LLFun,
                                TargetMean1,
                                TargetMean2,
                                TargetStd1,
                                TargetStd2,
                                Rho)


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
