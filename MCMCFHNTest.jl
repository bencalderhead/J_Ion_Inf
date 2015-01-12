using MCMCTypes

# Function for updating parameters - Gibbs step 1
include("MCMCUpdateParameters.jl")

# Function for the main MCMC routine
include("MCMCRun.jl")



# Define settings for simulation
OutputID            = "TestOutputFile"
Sampler             = "MH"
NumOfProposals      = 1
NumOfIterations     = 100
InitialStepSize     = 1.0
InitialiseFromPrior = false # Sample starting parameters from prior

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction = UpdateParameters


#######################
# Create an ODE model #
#######################

# Define the ODEModel object
ModelType    = "ODE"
ModelName    = "FithughNagumo"
NumOfParas   = 3
ParaNames    = ["a", "b", "c"]
DefaultParas = [0.2, 0.2, 3.0]

FHNPrior    = Array(PriorDistribution, 3)
FHNPrior[1] = PriorDistribution("Uniform", [0, 10])
FHNPrior[2] = PriorDistribution("Uniform", [0, 10])
FHNPrior[3] = PriorDistribution("Uniform", [0, 10])

# Data
DataMeasuments    = [1.0 1.0; 1.0 1.0; 1.0 1.0]
DataTimePoints    = [0.5 0.8; 0.5 0.8; 0.5 0.8]
DataNoiseVariance = [0.1 0.1; 0.1 0.1; 0.1 0.1]

# ODE specific
function FHN(t, y, ydot)
    a = 0.2
    b = 0.2
    c = 3.0
    ydot[1] = c*(y[1]-(y[1]^3)/3+y[2])
    ydot[2] = -(y[1]-a+b*y[2])/c
end

FHNFun                   = FHN
DefaultInitialConditions = [-1.0 1.0]
NumOfSpecies             = 3
ObservedSpecies          = [1 2 3]
UnobservedSpecies        = [0]
abstol                   = 1e-8
reltol                   = 1e-8
InferInitialConditions   = false


FHNModel = ODEModel( ModelType,
                     ModelName,
                     NumOfParas,
                     ParaNames,
                     DefaultParas,
                     UsePrior,
                     FHNPrior,
                     DataMeasuments,
                     DataTimePoints,
                     DataNoiseVariance,
                     FHNFun,
                     DefaultInitialConditions,
                     NumOfSpecies,
                     ObservedSpecies,
                     UnobservedSpecies,
                     abstol,
                     reltol,
                     InferInitialConditions)



############################
# Create simulation object #
############################

MySimulation = MCMCSimulation( OutputID,
                               Sampler,
                               NumOfProposals,
                               NumOfIterations,
                               InitialStepSize,
                               InitialiseFromPrior,
                               UpdateParasFunction,
                               FHNModel )


# Run the MCMC code, passing in the MCMCSimulation object
#RunMCMC( MySimulation )
