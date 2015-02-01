projbasedir = normpath(joinpath(dirname(Base.source_path()),"../.."))

#import needed modules
using Distributions
using Base.Test
using FactCheck

# Enum of Ion-channel models which can be created
AvailableModels = include(string(projbasedir,"ion/model/ModelEnum.jl"))

# Function to create available models
include(string(projbasedir,"ion/model/CreateModel.jl"))

push!(LOAD_PATH,projbasedir)
#need to add the project root to the path temperarily to import local MCMC types
using MCMCTypes
pop!(LOAD_PATH)

function TestIonModels()
    data = setup()
    TestBall(data)
    TestdCK(data)
end

function setup()
    #setup data
    data = Dict()
    data["burst_intervals"] = [0.1 0.2 0.1 0.2 0.15 0.16 0.18 0.05 0.1]
    data["burst_lengths"] = [3 1 5]
    data["tres"] = 2.5e-5
    data["tcrit"] = 1e-4
    data["useChs"] = 1
    data["conc"] = 1e-6
    return data
end

function TestdCK(data::Dict)
    LL = Array(Cdouble, 1)

    model = CreateModel(AvailableModels.DCK)

    testQ = model.GenerateQ([100.0; 100.0; 1000.0; 1000.0] , 1.0)
    realQ = [-100.0 100.0 0.0; 100.0 -1100.0 1000.0; 0.0 1000.0 -1000.0]

    Q = model.GenerateQ([100.0 100.0 1000.0 1e8], data["conc"])
    x = ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),LL, data["burst_intervals"] , length(data["burst_intervals"]),data["burst_lengths"],length(data["burst_lengths"]), Q , model.nopen, model.k, data["tres"], data["tcrit"], data["useChs"])

    facts("Checking del-Castillo-Katz") do
       
       context("Checking model fields") do
           @fact model.k => 3
           @fact model.DefaultParas => [100.0; 100.0; 1000.0; 1000.0]
           @fact model.UsePrior => true
           @fact model.NumOfParas => 4
           @fact model.nopen => 1
           @fact model.ParaNames => ["alpha1" "beta1" "km1" "kp1"]
       end

       context("Checking model priors") do
           AssociationRate = [false false false true] #which parameters are association rate constants
           for i=1:model.k
               @fact typeof(model.Priors[i]) => Uniform
               @fact model.Priors[i].a => 0.0
               if AssociationRate[i]
                   @fact model.Priors[i].b => 1e10
               else
                   @fact model.Priors[i].b => 1e6
               end
           end
       end

       context("Checking Q matrix") do
           for i=1:model.k
               for j=1:model.k
                   @fact testQ[i,j] => roughly(realQ[i,j]; atol = 1e-10) "Mismatch at $i,$j"
               end
           end
       end
       
       context("Checking likelihood") do
           @fact LL[1] => roughly(-52.47666128726489; atol =1e-10 )
       end
    end
end

function TestBall(data::Dict)
    LL = Array(Cdouble, 1)

    #Test 1: Ball Ion model
    model = CreateModel(AvailableModels.BALL1989)

    testQ = model.GenerateQ([100.0; 5e6], 1.0)
    realQ = [ -100.0 100.0; 5e6 -5e6]

    Q = model.GenerateQ([100.0; 5e6], data["conc"])
    x = ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),LL, data["burst_intervals"] , length(data["burst_intervals"]),data["burst_lengths"],length(data["burst_lengths"]), Q , model.nopen, model.k, data["tres"], data["tcrit"], data["useChs"])

    facts("Checking Ball 1989") do
       
       context("Checking model fields") do
           @fact model.k => 2
           @fact model.DefaultParas => [100.0; 5e6]
           @fact model.UsePrior => true
           @fact model.NumOfParas => 2
           @fact model.nopen => 1
           @fact model.ParaNames => ["alpha1" "beta1"]
       end

       context ("Checking first prior") do
           @fact typeof(model.Priors[1]) => Uniform
           @fact model.Priors[1].a => 0.0
           @fact model.Priors[1].b => 1e6
       end
       
       context("Checking model priors") do
           AssociationRate = [false, true] #which parameters are association rate constants
           for i=1:model.k
               @fact typeof(model.Priors[i]) => Uniform
               @fact model.Priors[i].a => 0.0
               if AssociationRate[i]
                   @fact model.Priors[i].b => 1e10
               else
                   @fact model.Priors[i].b => 1e6
               end
           end
       end

       context("Checking Q matrix") do
           for i=1:model.k
               for j=1:model.k
                   @fact testQ[i,j] => roughly(realQ[i,j]; atol = 1e-10) "Mismatch at $i,$j"
               end
           end
       end
       
       context("Checking likelihood") do
           @fact LL[1] => -52.56916031035307
       end

    end
end

