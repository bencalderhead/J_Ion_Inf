using Distributions
using MCMCTypes
using Sundials

# Function for updating parameters - Gibbs step 1
include("MCMCUpdateParameters.jl")

# Function for updating indicator - Gibbs step 2
include("MCMCUpdateIndicator.jl")

# Function for the main MCMC routine
include("MCMCRun.jl")



# Define settings for simulation
#Sampler             = "MH"
#Sampler             = "AdaptiveMH"
#Sampler             = "SmMALA"
Sampler             = "TrSmMALA"

if Sampler == "MH"
  OutputID            = "./Output/PredatorPrey_MH_Output"
  InitialStepSize     = 0.02 # MH stepsize
  NumOfIterations     = 15000
elseif Sampler == "AdaptiveMH"
  OutputID            = "./Output/PredatorPrey_AdaptiveMH_Output"
  InitialStepSize     = 0.0001 # MH stepsize
  NumOfIterations     = 30000
elseif Sampler == "SmMALA"
  OutputID            = "./Output/PredatorPrey_SmMALA_Output"
  InitialStepSize     = 0.7 # SmMALA stepsize
  NumOfIterations     = 15000
elseif Sampler == "TrSmMALA"
  OutputID            = "./Output/PredatorPrey_TrSmMALA_Output"
  InitialStepSize     = 100 # TrSmMALA stepsize
  NumOfIterations     = 5000
end

NumOfProposals      = 1
ProposalCovariance  = eye(6)
InitialiseFromPrior = false # Sample starting parameters from prior

AuxiliaryVars       = 20; # Number of tangent samples (We'll store this in the model for now)

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction     = UpdateParameters

# Define function for sampling the indicator variable (2nd Gibbs step)
SampleIndicatorFunction = SampleIndicator


#######################
# Create an ODE model #
#######################

# Define the Model object
ModelType    = "ODE"
ModelName    = "PredatorPrey"
NumOfParas   = 6
ParaNames    = ["r" "k" "s" "a" "u" "v"]
#DefaultParas = [0.6; 100; 1.2; 25; 0.5; 0.3] # Column vector
DefaultParas = [0.358755779146462; 107.376260452987; 0.872321393358852; 53.3621151904749; 0.637996153522885; 0.272318515919793]

UsePrior     = true
Prior        = Array(Distribution, NumOfParas)
for i = 1:NumOfParas
    Prior[i] = Uniform(0, 150)
end

# Specify the data
DataMeasurements = [51.7002527839933	3.62886163423551;
88.4059566464213	11.8592208786134;
50.5971643551162	33.1299499686354;
25.2455746232269	30.1443180597855;
42.2660932030265	9.00691042344652;
51.4088102060744	29.8816967299902]

DataTimePoints = [linspace(0,50,6) linspace(0,50,6)]


DataNoiseVariance = ones(6,2) # Same size as the dataset
DataNoiseVariance[:,1] = sqrt(10)
DataNoiseVariance[:,2] = sqrt(10)

NumOfTimePoints = 6

# Specify differential equation function
function PredatorPrey(t,y,ydot,p)
    P = y[1]
    Q = y[2]

    r = p[1]
    K = p[2]
    s = p[3]
    a = p[4]
    u = p[5]
    v = p[6]

    tmp = P*Q/(a + P);

    ydot[1] = r*P*(1-P/K) - s*tmp;
    ydot[2] = u*tmp - v*Q;

    nothing
end


ODEFunction              = PredatorPrey
DefaultInitialConditions = [50.0, 5.0]
NumOfSpecies             = 2
ObservedSpecies          = [1, 2]
UnobservedSpecies        = []
abstol                   = 1e-6
reltol                   = 1e-6
InferInitialConditions   = false

MyODEModel = ODEModel( ModelType,
                       ModelName,
                       NumOfParas,
                       ParaNames,
                       DefaultParas,
                       UsePrior,
                       Prior,
                       DataMeasurements,
                       DataTimePoints,
                       DataNoiseVariance,
                       NumOfTimePoints,
                       ODEFunction,
                       DefaultInitialConditions,
                       NumOfSpecies,
                       ObservedSpecies,
                       UnobservedSpecies,
                       abstol,
                       reltol,
                       InferInitialConditions,
                       AuxiliaryVars)


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
                               MyODEModel )


# Run the MCMC code, passing in the MCMCSimulation object
MCMCRun( MySimulation )
