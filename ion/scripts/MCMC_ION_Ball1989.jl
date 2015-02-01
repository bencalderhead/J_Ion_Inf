
#we need to key off the absolute path of the project root so to include the relevent files
projbasedir = normpath(joinpath(dirname(Base.source_path()),"../.."))

#import needed modules
using Distributions
using MAT

push!(LOAD_PATH,projbasedir)
#need to add the project root to the path temperarily to import local MCMC types
using MCMCTypes
pop!(LOAD_PATH)

# Function for updating parameters - Gibbs step 1
include(string(projbasedir,"MCMCUpdateParameters.jl"))

# Function for updating indicator - Gibbs step 2
include(string(projbasedir,"MCMCUpdateIndicator.jl"))

# Function for the main MCMC routine
include(string(projbasedir,"MCMCRun.jl"))

# Enum of Ion-channel models which can be created
AvailableModels = include(string(projbasedir,"ion/model/ModelEnum.jl"))

# Function to create available models
include(string(projbasedir,"ion/model/CreateModel.jl"))

# Function to parse matlab data files
include(string(projbasedir,"ion/data/MatlabParser.jl"))


# Function to create an ion-channel model

# Define settings for simulation, only MH implemented for models of Type ION
#Sampler             = "MH"
Sampler             = "AdaptiveMH"
#Sampler             = "SmMALA"
#Sampler             = "TrSmMALA"

if Sampler == "MH"
  OutputID            = "$projbasedir/Output/ION_Ball_MH_Output"
  InitialStepSize     = 1 # MH stepsize
  NumOfIterations     = 15000
elseif Sampler == "AdaptiveMH"
  OutputID            = "$projbasedir/Output/ION_Ball_AdaptiveMH_Output"
  InitialStepSize     = 1 # MH stepsize
  NumOfIterations     = 10000
elseif Sampler == "SmMALA"
  OutputID            = "$projbasedir/Output/ION_Ball_SmMALA_Output"
  InitialStepSize     = 0.7 # SmMALA stepsize
  NumOfIterations     = 15000
elseif Sampler == "TrSmMALA"
  OutputID            = "$projbasedir/Output/ION_Ball_TrSmMALA_Output"
  InitialStepSize     = 100 # TrSmMALA stepsize
  NumOfIterations     = 5000
end

NumOfProposals      = 1
ProposalCovariance  = [ 10 0; 0 1000000 ]
InitialiseFromPrior = false # Sample starting parameters from prior

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction     = UpdateParameters

# Define function for sampling the indicator variable (2nd Gibbs step)
SampleIndicatorFunction = SampleIndicator

#create the Ball model
println("Creating model Ball1989")
model = CreateModel(AvailableModels.BALL1989)

#data parsing from MATLAB experimental files

datafile = string(projbasedir,"ion/data/balldata.mat")
println("Parsing datafile $datafile")
data = MatlabParser(datafile)

println(string("Likelihood with default params = ", model.LLEval(model.DefaultParas,model,data)))

println(string("Likelihood with MAP params = ", model.LLEval([1003.4028; 9962666.0321],model,data)))

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
                               model,
                               data )


# Run the MCMC code, passing in the MCMCSimulation object
MCMCRun( MySimulation )

