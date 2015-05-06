module MCMCTypes

export IONModel, ODEModel, TargetOnlyModel, GaussianBivariate, MarkovChain, MarkovChainGeometry, MarkovChainProposal, MCMCSimulation
export ProposalDistributionMH, ProposalDistributionSmMALA, ProposalDistributionSmMALARandom, ProposalDistributionAdaptiveMH

# Define the statistical model object
type IONModel
    #fields from Target Only Model, acts as interface to MCMCRun.jl

    ModelType
    ModelName
    NumOfParas::Int64
    ParaNames::Array
    DefaultParas::Array{Float64}
    UsePrior::Bool
    Priors::Array{Any} #types from the distribution Pkg
    LLEval::Function

    #Data
    #NumberOfConcs::Int64
    #Concs::Array{Float64} #concentrations at each recording
    #Burstnos::Array{Int64} # number of bursts in each concentration
    #BurstLengths::Array{Any} #numberofconcs in length, each index holds an Array{Int64} of the number of intervals in each burst
    #BurstIntervals::Array{Any} #numberofconcs in length, each index holds an Array{Float64} which are the intervals across all bursts. burstlengths is used to delineate each burst in dcprogs_wrapper
    #Tres::Array{Float64} #resolution time for each set
    #Tcrit::Array{Float64} #time interval to delineate bursts
    #useChs::Array{Bool}   #whether to use CHS vectors for initial occupacies

    #ION specific
    GenerateQ::Function #generates a Q-matrix which encodes the transition rates of the continuous time Markov process
    nopen::Int64
    k::Int64

end

# Define the statistical model object
type ODEModel

    ModelType
    ModelName
    NumOfParas::Int64
    ParaNames::Array
    DefaultParas::Array{Float64}
    UsePrior::Bool
    Priors::Array{Any}

    # Data
    DataMeasurements::Array{Float64}
    DataTimePoints::Array{Float64}
    DataNoiseVariance::Array{Float64}
    NumOfTimePoints::Int64

    # ODE specific
    ODEFunction::Function
    DefaultInitialConditions::Array{Float64}
    NumOfSpecies::Int64
    ObservedSpecies::Array{Int64}
    UnobservedSpecies::Array{Int64}
    abstol::Float64
    reltol::Float64
    InferInitialConditions::Bool

    AuxiliaryVars::Any

end


# Define the statistical model object
type GaussianBivariate

    ModelType
    ModelName
    NumOfParas::Int64
    ParaNames::Array
    DefaultParas::Array{Float64}
    UsePrior::Bool
    Priors::Array{Any}
    LLEval::Function

    TargetMean1::Float64
    TargetMean2::Float64
    TargetStd1::Float64
    TargetStd2::Float64
    Rho::Float64
    TargetCov::Array{Float64}

    # Now set up the functions
    function GaussianBivariate(ModelType, ModelName, NumOfParas, ParaNames, DefaultParas, UsePrior, Priors, LLEval, TargetMean1, TargetMean2, TargetStd1, TargetStd2, Rho)
        this              = new()
        this.ModelType    = ModelType
        this.ModelName    = ModelName
        this.NumOfParas   = NumOfParas
        this.ParaNames    = ParaNames
        this.DefaultParas = DefaultParas
        this.UsePrior     = UsePrior
        this.Priors       = Priors
        this.LLEval       = LLEval

        this.TargetMean1  = TargetMean1
        this.TargetMean2  = TargetMean2
        this.TargetStd1   = TargetStd1
        this.TargetStd2   = TargetStd2
        this.Rho          = Rho
        this.TargetCov    = [TargetStd1^2 TargetStd1*TargetStd2*Rho; TargetStd1*TargetStd2*Rho TargetStd2^2]

        return this
    end
end


# Define the statistical model object
type TargetOnlyModel

    ModelType
    ModelName
    NumOfParas::Int64
    ParaNames::Array
    DefaultParas::Array{Float64}
    UsePrior::Bool
    Priors::Array{Any}
    LLEval::Function
end


# Define Markov chain geometry object
type MarkovChainGeometry

    Parameters::Array{Float64}

    LL::Float64
    LogPrior::Float64
    ProposalProbability::Array{Float64}

    GradLL::Array{Float64}              # vector
    GradLogPrior::Array{Float64}        # vector
    HessianLL::Array{Float64}           # matrix
    HessianGradLogPrior::Array{Float64} # matrix

end


# Define proposal object for M-H MCMC
type ProposalDistributionMH

    Density::Any                        # This will be a distribution type
    StepSize::Float64
    ProposalCov::Array{Float64}

end

# Define proposal object for SmMALA MCMC
type ProposalDistributionSmMALA

    Density::Any                        # This will be a distribution type
    StepSize::Float64

end

# Define proposal object for SmMALA MCMC
type ProposalDistributionSmMALARandom

    Density::Any                        # This will be a distribution type
    StepSize::Float64
    NumberOfVectors::Int64
    TangentVectors::Array{Float64}

end

# Define proposal object for Adaptive MCMC
type ProposalDistributionAdaptiveMH

    Density::Any                        # This will be a distribution type
    StepSize::Float64
    RunningMean::Array{Float64}
    RunningCov::Array{Float64}

end


# Define Markov chain object
type MarkovChain

    Sampler
    NumOfProposals::Int64
    SampleIndicator::Int64
    Initialised::Bool
    CurrentIteration::Int64

    AttemptedProposal::Int64
    AcceptedProposal::Int64
    AttemptedExchange::Int64
    AcceptedExchange::Int64

    Geometry::Array{MarkovChainGeometry}

    ProposalDistribution::Any           # Distribution for the proposal
end


# Define the MCMC Simulation object
type MCMCSimulation

    OutputID
    Sampler # e.g. MH, MALA, SmMALA, TrSmMALA etc
    NumOfProposals::Int64
    NumOfIterations::Int64
    InitialStepSize::Float64
    ProposalCovariance::Array{Float64}
    InitialiseFromPrior::Bool

    # Function for updating parameter values from proposal kernel
    UpdateParameters::Function

    # Function for sampling the indicator variable
    SampleIndicator::Function

    # This can be any model, as long as corresponding UpdateParameters function is defined
    Model
    
    # This is a Dictionary which the LLEval function in the model needs to understand
    Data
end


end
