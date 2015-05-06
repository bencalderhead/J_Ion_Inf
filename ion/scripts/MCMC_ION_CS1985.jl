@everywhere projbasedir = "/home/michaelepstein/Ion/"

#import needed modules
@everywhere using Distributions
@everywhere using MAT

@everywhere push!(LOAD_PATH,projbasedir)
#need to add the project root to the path temperarily to import local MCMC types
@everywhere using MCMCTypes
@everywhere pop!(LOAD_PATH)

# Function for updating parameters - Gibbs step 1
@everywhere include(string(projbasedir,"MCMCUpdateParameters.jl"))

# Function for updating indicator - Gibbs step 2
@everywhere include(string(projbasedir,"MCMCUpdateIndicator.jl"))

# Function for the main MCMC routine
@everywhere include(string(projbasedir,"MCMCRun.jl"))

# Enum of Ion-channel models which can be created
@everywhere AvailableModels = include(string(projbasedir,"ion/model/ModelEnum.jl"))

# Function to create available models
@everywhere include(string(projbasedir,"ion/model/CreateModel.jl"))

# Function to parse matlab data files
@everywhere include(string(projbasedir,"ion/data/MatlabParser.jl"))

@everywhere f(p) = model.LLEval(p,model,data)




# Define settings for simulation, only MH implemented for models of Type ION
Sampler             = "MH"
#Sampler             = "AdaptiveMH"

if Sampler == "MH"
  OutputID            = "./Output/ION_CS1985_MH_Output"
  InitialStepSize     = 1000 # MH stepsize
  NumOfIterations     = 15
elseif Sampler == "AdaptiveMH"
  OutputID            = "./Output/ION_CS1985_AdaptiveMH_Output"
  InitialStepSize     = 0.0001 # MH stepsize
  NumOfIterations     = 15
end

NumOfProposals      = 1000
ProposalCovariance  = eye(9)
InitialiseFromPrior = false # Sample starting parameters from prior

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction     = UpdateParameters

# Define function for sampling the indicator variable (2nd Gibbs step)
SampleIndicatorFunction = SampleIndicator


#data parsing from MATLAB experimental files
datafile = "/home/michaelepstein/Utsire/data/AchRealData.mat"
println("Parsing datafile $datafile")
@everywhere data = MatlabParser(datafile)

#create the ion model
@everywhere model = CreateModel(AvailableModels.CS1985)
println(string("Likelihood with default params = ", model.LLEval(DefaultParas,model)))

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
                               model, )
                               data

# Run the MCMC code, passing in the MCMCSimulation object
MCMCRun( MySimulation )

