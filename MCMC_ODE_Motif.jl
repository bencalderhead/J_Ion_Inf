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
Sampler             = "SmMALA"
#Sampler             = "TrSmMALA"

if Sampler == "MH"
  OutputID            = "./Output/Motif_MH_Output_Start2"
  InitialStepSize     = 0.02 # MH stepsize
  NumOfIterations     = 15000
elseif Sampler == "AdaptiveMH"
  OutputID            = "./Output/Motif_AdaptiveMH_Output_Start1"
  InitialStepSize     = 0.0001 # MH stepsize
  NumOfIterations     = 15000
elseif Sampler == "SmMALA"
  OutputID            = "./Output/Motif_SmMALA_Output"
  InitialStepSize     = 0.0001 # SmMALA stepsize
  NumOfIterations     = 15000
elseif Sampler == "TrSmMALA"
  OutputID            = "./Output/Motif_TrSmMALA_Output_Start2"
  InitialStepSize     = 0.5 # TrSmMALA stepsize
  NumOfIterations     = 5000
end

NumOfProposals      = 1
ProposalCovariance  = eye(7)
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
ModelName    = "Motif"
NumOfParas   = 7
ParaNames    = ["k1" "k2" "k3" "k4" "k5" "k6" "k7"]
#DefaultParas = [1.0; 1.0; 1.0; 1.0; 1.0; 1.0; 1.0] # Column vector
DefaultParas = [1.0; 0.9; 0.7; 1.0; 0.6; 5.6; 8.5] # Column vector
#DefaultParas = [1.0; 0.9; 0.7; 1.0; 9.7; 5.6; 0.6] # Column vector

UsePrior     = true
Prior        = Array(Distribution, NumOfParas)
for i = 1:NumOfParas
    Prior[i] = Uniform(0, 10)
end

# Specify the data
DataMeasurements = [1	0	1	0	0	1	0;
0.445196174826338	0.346866344955917	0.773107204113088	0.163440412691046	0.0484444535780175	0.991608497937154	0.00800432304644002;
0.239580636452932	0.508512359351390	0.745204298222111	0.134103878318130	0.0979802690198388	0.972557634753506	0.0281674818119963;
0.0630032656157987	0.623817267742466	0.764094992376031	0.0584640757989290	0.120325571268486	0.945954628816380	0.0547306006925542;
0.0273147186782370	0.685082981942735	0.821906992902959	0.0256695597435487	0.109665901541326	0.945459571620418	0.0533775264499349;
0.00694543265934263	0.685514112180426	0.867255042040815	0.0136238326067817	0.0784932647571628	0.954943315938326	0.0437745479285641;
0.0150012916598958	0.693056415626388	0.911477804097201	0.00345762326160713	0.0591629494513626	0.967609322258256	0.0345641660437732]

DataTimePoints = [0 0.5 1 2 3 4 5;
                  0 0.5 1 2 3 4 5;
                  0 0.5 1 2 3 4 5;
                  0 0.5 1 2 3 4 5;
                  0 0.5 1 2 3 4 5;
                  0 0.5 1 2 3 4 5;
                  0 0.5 1 2 3 4 5]

DataTimePoints = DataTimePoints'

DataNoiseVariance = ones(7,7) # Same size as the dataset
DataNoiseVariance[:,1] = 0.000331408785760581
DataNoiseVariance[:,2] = 0.000172298359044859
DataNoiseVariance[:,3] = 2.13620891473913e-05
DataNoiseVariance[:,4] = 1.15808937452288e-05
DataNoiseVariance[:,5] = 4.26512029018256e-06
DataNoiseVariance[:,6] = 1.09623757850807e-06
DataNoiseVariance[:,7] = 1.09623757850806e-06

NumOfTimePoints = 7

# Specify differential equation function
function Motif(t,y,ydot,p)
    S      = y[1]
    dS     = y[2]
    R      = y[3]
    RS     = y[4]
    Rpp    = y[5]
    PhA    = y[6]
    RppPhA = y[7]

    k1 = p[1]
    k2 = p[2]
    k3 = p[3]
    k4 = p[4]
    k5 = p[5]
    k6 = p[6]
    k7 = p[7]

    dS      = -k1*S - k2*S*R + k3*RS
    ddS     = k1*S
    dR      = -k2*S*R + k3*RS + k7*RppPhA
    dRS     = k2*S*R - k3*RS - k4*RS
    dRpp    = k4*RS - k5*Rpp*PhA + k6*RppPhA
    dPhA    = -k5*Rpp*PhA + k6*RppPhA + k7*RppPhA
    dRppPhA = k5*Rpp*PhA - k6*RppPhA - k7*RppPhA

    ydot[1] = dS
    ydot[2] = ddS
    ydot[3] = dR
    ydot[4] = dRS
    ydot[5] = dRpp
    ydot[6] = dPhA
    ydot[7] = dRppPhA
end


ODEFunction              = Motif
DefaultInitialConditions = [1., 0., 1., 0., 0., 1., 0.]
NumOfSpecies             = 7
ObservedSpecies          = [1, 2, 3, 4]
UnobservedSpecies        = [5, 6, 7]
abstol                   = 1e-8
reltol                   = 1e-8
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
