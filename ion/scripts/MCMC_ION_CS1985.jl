using Distributions
using MCMCTypes
using MAT

# Function for updating parameters - Gibbs step 1
include("MCMCUpdateParameters.jl")

# Function for updating indicator - Gibbs step 2
include("MCMCUpdateIndicator.jl")

# Function for the main MCMC routine
include("MCMCRun.jl")

# Define settings for simulation, only MH implemented for models of Type ION
Sampler             = "AdaptiveMH"
#Sampler             = "AdaptiveMH"
#Sampler             = "SmMALA"
#Sampler             = "TrSmMALA"

if Sampler == "MH"
  OutputID            = "./Output/ION_CS1985_MH_Output"
  InitialStepSize     = 1000 # MH stepsize
  NumOfIterations     = 15000
elseif Sampler == "AdaptiveMH"
  OutputID            = "./Output/ION_CS1985_AdaptiveMH_Output"
  InitialStepSize     = 0.0001 # MH stepsize
  NumOfIterations     = 30000
elseif Sampler == "SmMALA"
  OutputID            = "./Output/ION_CS1985_SmMALA_Output"
  InitialStepSize     = 0.7 # SmMALA stepsize
  NumOfIterations     = 15000
elseif Sampler == "TrSmMALA"
  OutputID            = "./Output/ION_CS1985_TrSmMALA_Output"
  InitialStepSize     = 100 # TrSmMALA stepsize
  NumOfIterations     = 5000
end

NumOfProposals      = 1
ProposalCovariance  = eye(9)
InitialiseFromPrior = false # Sample starting parameters from prior

# Define function for updating the parameters (1st Gibbs step)
UpdateParasFunction     = UpdateParameters

# Define function for sampling the indicator variable (2nd Gibbs step)
SampleIndicatorFunction = SampleIndicator

#define the model

ModelType    = "ION Acetylcholine channel model"
ModelName    = "C&S 1985"
NumOfParas   = 9
ParaNames    = ["kp2*" "alpha1" "alpha2" "beta2" "km2" "beta1" "kp2" "km1" "kp1"] #km2* set by mr
DefaultParas = [50/1e-7; 3000; 500; 15000; 2000; 15; 50/1e-7; 2000; 5/1e-7]

UsePrior     = true
Prior        = Array(Distribution, NumOfParas)
for i = 1:NumOfParas
    Prior[i] = Uniform(0, 10e10)
end

#ion channel specific params and functions
nopen = 2
k = 5

function generateQ(params::Array{Float64},conc::Float64)
  kp2s = params[1]
  α1 = params[2]
  α2 = params[3]
  β2 = params[4]
  km2 = params[5]
  β1 = params[6]
  kp2 = params[7]
  km1 = params[8]
  kp1 = params[9]

  km2s = (kp2s * α2 * km2 * β1 )/(α1 * kp2 * β2); #parameter set by mr

  Q = [ [-(α1+(conc*kp2s))  (conc*kp2s)         0                α1                         0;]   
        [2*km2s             -(α2 + (2 * km2s))  α2               0                          0;]
        [0                  β2                  -(β2 + (2*km2))  2*km2                      0;]
        [β1                 0                   kp2*conc         -(β1 + (kp2*conc) + km1)   km1;]
        [0                  0                   0               2*kp1*conc                  -2*kp1*conc;]
      ]

end

function LLFunction(params::Array{Float64}, model::IONModel)

  LL = 0 
  SetLikelihood = Array(Cdouble,1)
  for i=1:model.NumberOfConcs
    Q = model.GenerateQ(params,model.Concs[i])
 
    x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),SetLikelihood,model.BurstIntervals[i],length(model.BurstIntervals[i]),model.BurstLengths[i],length(model.BurstLengths[i]), Q , model.nopen, model.k, model.Tres[i],model.Tcrit[i],model.useChs[i])
    LL += SetLikelihood[1]
  end
  return LL
end

#data parsing from MATLAB experimental files

datafile = "/home/michaelepstein/Utsire/data/AchRealData.mat"
file = matopen(datafile)

bursts = read(file,"bursts")
concs = read(file,"concs")
tcrit = read(file,"tcrit")
tres = read(file,"tres")
useChs =  read(file,"useChs")
close(file)

experimentno = length(concs)
println("number of sets = $experimentno")
burstnos = Array(Csize_t,experimentno)
bsts = Array(Any,experimentno)
burst_lengths = Array(Any,experimentno)

#need the indexes for each burst set to delineate bursts
counter=1
for i=1:experimentno
  burstnos[i] = length(bursts[i])
  burst_lengths[i] = zeros(Csize_t,burstnos[i])
  for j=1:length(bursts[i])
    burst_lengths[i][j] = length(bursts[i][j])
  end

  bsts[i] = zeros(Float64,sum(burst_lengths[i]))
  counter = 1
  for j = 1:length(bursts[i])
    for m = 1:length(bursts[i][j])
      bsts[i][counter] = bursts[i][j][m]
      counter += 1
    end
  end
  println(string("number of burst intervals = " ,length(bsts[i])))
  println(string("number of bursts = ", sum(burstnos[i])))
end

#create the ion model
model = IONModel(ModelType,ModelName,NumOfParas,ParaNames,DefaultParas,UsePrior,Prior,LLFunction,experimentno,concs,burstnos,burst_lengths,bsts,tres,tcrit,useChs,generateQ,nopen,k)

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
                               model )


# Run the MCMC code, passing in the MCMCSimulation object
MCMCRun( MySimulation )




