#include the enum based on the absolute path of the source file
using PyCall
using Distributions
@pyimport dcprogs.likelihood as dcplikelihood
modeltypes = include(string(dirname(Base.source_path()),"/ModelEnum.jl"))
function CreateModel(x::Int64,useC=true)
    #function creates a model as defined in ModelEnum.jl

    if (x == modeltypes.BALL1989)
        println("Creating Ball 1989 Model. Params = 2, States = 2")
        ModelType    = "ION TwoState channel model"
        ModelName    = "Ball 1989"
        NumOfParas   = 2
        ParaNames    = ["alpha1" "beta1"] #km2* set by mr
        DefaultParas = [100.0 ; 5e6]
        ConcDependent= [false, true] #local variable used to set priors on rate constants
        UsePrior     = true
        Prior        = Array(Distribution, NumOfParas)

        for i = 1:NumOfParas
            if ConcDependent[i]
                Prior[i] = Uniform(0, 1e10)
            else
                Prior[i] = Uniform(0, 1e6)
            end
        end
        #ion channel specific params and functions
        nopen = 1
        k = 2

        function generateQ(params::Array{Float64},conc::Float64)
            α1 = params[1]
            β1 = params[2]

            Q = [ [-α1     α1;        ]
                [β1*conc -(β1*conc); ]
            ]
        end


    elseif (x == modeltypes.DCK)
        println("Creating del-Castillo-Katz Model. Params = 4, States = 3")
        ModelType    = "ION Three State channel model"
        ModelName    = "del-Castillo Katz"
        NumOfParas   = 4
        ParaNames    = ["alpha1" "beta1" "km1" "kp1"]
        ConcDependent= [false,false,false,true]
        DefaultParas = [100.0 ; 100.0; 1000.0; 1000.0]

        UsePrior     = true
        Prior        = Array(Distribution, NumOfParas)
        for i = 1:NumOfParas
            if ConcDependent[i]
                Prior[i] = Uniform(0, 1e10)
            else
                Prior[i] = Uniform(0, 1e6)
            end
        end

        #ion channel specific params and functions
        nopen = 1
        k = 3

       function generateQ(params::Array{Float64},conc::Float64)
           α1 = params[1]
           β1 = params[2]
           km1 = params[3]
           kp1 = params[4]

           Q = [ [-α1    α1         0;]
               [β1     -(β1+km1)  km1; ]
               [0      kp1*conc   -kp1*conc;]
           ]
       end


    elseif (x == modeltypes.CS1985)
        println("Creating Colquhoun and Sakmann 1985 Model. Params = 9, States = 5")
        error("Model not implemented yet - no model created")

    elseif (x == modeltypes.HATTON2003)
        println("Creating Hatton 2003 Model with independent binding sites. Params = 10, States = 7")
        error("Model not implemented yet - no model created")

    else
        error("Unknown model type for creation - no model created")
    end

    #likelihood function identical for all models
    function LLFunction(params::Array{Float64}, model::IONModel, data::Dict)
      
      LL = 0
      SetLikelihood = Array(Cdouble,1)
      for i=1:data["experiment_nos"]
            Q = model.GenerateQ(params,data["concs"][i])

            x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Csize_t,Ptr{Csize_t},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),SetLikelihood,data["burst_intervals"][i],length(data["burst_intervals"][i]),data["burst_lengths"][i],length(data["burst_lengths"][i]), Q , model.nopen, model.k, data["tres"][i],data["tcrit"][i],data["useChs"][i])
            LL += SetLikelihood[1]
        end
        return LL
    end

    function pyLLFunction(params::Array{Float64}, model::IONModel, data::Dict)
      LL = 0

      for i=1:data["experiment_nos"]
        Q = model.GenerateQ(params,data["concs"][i])
        pylik = dcplikelihood.Log10Likelihood(data["pybursts"][i], model.nopen, data["tres"][i], data["tcrit"][i])
        pyQmat = dcplikelihood.QMatrix(Q , model.nopen)
        LL+=pycall(pylik,PyAny,pyQmat)*log(10)
      end
      return LL
    end

    #create the model here and return it. Choose to use the C wrapper or the python
    if useC
        model = IONModel(ModelType,ModelName,NumOfParas,ParaNames,DefaultParas,UsePrior,Prior,LLFunction,generateQ,nopen,k) 
    else

        model = IONModel(ModelType,ModelName,NumOfParas,ParaNames,DefaultParas,UsePrior,Prior,pyLLFunction,generateQ,nopen,k) 

    end


end
