using FactCheck
using PyCall
@pyimport dcprogs.likelihood as dcplikelihood

function setup()
    nopen = 2
    k = 5

    tres = 5e-5

    Q = [ [ -3050 50           0      3000   0; ]
          [ 2/3   -(500 + 2/3) 500    0      0; ]
          [ 0     15000        -19000 4000   0; ]
          [ 15    0            50     -2065  2000;]
          [ 0     0            0      10     -10;] ]

    #standard dcp options
    dcProgsOptions = [2 1e-12 1e-12 100 -1e6 0] 

    t = [1e-5 2e-5 3e-5 4e-5 5e-5 6e-5 7e-5 8e-5 9e-5 1e-4]

    data = Dict()
    data["Q"] = Q
    data["tres"] = tres
    data["nopen"] = nopen
    data["k"] = k
    data["dcProgsOptions"] = dcProgsOptions
    data["t"] = t
    data["tcrit"] = 1e-4

    return data
end

function TestIonPredictive()
    #setup code to match DCPROGS
    println("Running tests from CHS 1996...")
    data = setup()
    TestOccupancies(data)
    TestIdealAFt(data)
    TestIdealFAt(data)

    #TestIdealOpenTime(data)
    #TestIdealCloseTime(data) 
end

function TestIdealCloseTime(data::Dict)

    idealPdf = zeros(length(data["t"]))
    FAt = zeros(3 , 2)
    phiF = zeros (1 , 3)
    perfectres = 0.0

    x=ccall((:dcpOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool, Ptr{Float64}), phiF, data["Q"] , data["nopen"] , data["k"], perfectres, false, data["dcProgsOptions"] )

    for i=1:length(data["t"])
        x=ccall((:dcpIdealGXYt,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool ), FAt, data["Q"] , 2 , 5, data["t"][i], false )
        idealPdf[i] = (data["t"][i] * phiF * FAt * ones(data["nopen"]))[1]
    end

    facts("Ideal Open Time Pdf") do
        @fact size(idealPdf) => (10,)
        @pending idealPdf => roughly([0.0067623400463521785,0.013349872529189111,0.019768740249407404,0.026024851516588458,0.03212388888198266,0.03807131755712209,0.04387239352907677,0.0495321713829914,0.05505551184217408,0.06044708903565556],atol=1e-12)
    end
end

function TestIdealOpenTime(data::Dict)

    idealPdf = zeros(length(data["t"]))
    AFt = zeros(2,3)
    phiA = zeros(1,2)
    perfectres = 0.0

    x=ccall((:dcpOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool, Ptr{Float64}), phiA, data["Q"] , data["nopen"] , data["k"], perfectres, true, data["dcProgsOptions"] )

    for i=1:length(data["t"])
        x=ccall((:dcpIdealGXYt,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool ), AFt, data["Q"] , 2 , 5, data["t"][i], true )
        idealPdf[i] = (data["t"][i] * phiA * AFt * ones(data["k"]-data["nopen"]))[1]
    end

    facts("Ideal Open Time Pdf") do
        @fact size(idealPdf) => (10,)
        @pending idealPdf => roughly([0.0067623400463521785,0.013349872529189111,0.019768740249407404,0.026024851516588458,0.03212388888198266,0.03807131755712209,0.04387239352907677,0.0495321713829914,0.05505551184217408,0.06044708903565556],atol=1e-12)
    end
end

function TestOccupancies(data::Dict)
   #dcpOccupancies ( double * phi, double * Q, size_t nopen, size_t k, double tau, bool initial, double * dcpOptions )
   phiA = zeros(1, data["nopen"])
   x=ccall((:dcpOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool, Ptr{Float64}), phiA, data["Q"] , data["nopen"] , data["k"], data["tres"], true, data["dcProgsOptions"] )
   
   phiF = zeros(data["k"] - data["nopen"])
   x=ccall((:dcpOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool, Ptr{Float64}), phiF, data["Q"] , data["nopen"] , data["k"], data["tres"], false, data["dcProgsOptions"] ) 

   chsPhiA = zeros(1, data["nopen"])
   x=ccall((:dcpCHSOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Float64,Bool, Ptr{Float64}), chsPhiA, data["Q"] , data["nopen"] , data["k"], data["tres"], data["tcrit"], true, data["dcProgsOptions"] )


   chsPhiF = zeros(data["k"] - data["nopen"])
   x=ccall((:dcpCHSOccupancies,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Float64,Bool, Ptr{Float64}), chsPhiF, data["Q"] , data["nopen"] , data["k"], data["tres"], data["tcrit"], false, data["dcProgsOptions"] )

   pyQmat = dcplikelihood.QMatrix(data["Q"] , data["nopen"])
   missedg = dcplikelihood.MissedEventsG(pyQmat,data["tres"],nmax=data["dcProgsOptions"][1],xtol=data["dcProgsOptions"][2],rtol=data["dcProgsOptions"][3],itermax=data["dcProgsOptions"][4],lower_bound=data["dcProgsOptions"][5],upper_bound=data["dcProgsOptions"][6])

   facts("Checking Missed Events Occupancies") do
       context("Checking Missed Events occupancies") do
           @fact size(phiA) => (1,2)
           @fact phiA => missedg[:initial_occupancies]
           @fact size(phiF) => (3,)
           @fact phiF => missedg[:final_occupancies]
       end
       context("Checking CHS occupances") do
           @fact size(chsPhiA) => (1,2)
           @fact chsPhiA => missedg[:initial_CHS_occupancies](data["tcrit"]) "dcprogs python gives "
           @fact size(chsPhiF) => (3,)
           @fact chsPhiF => missedg[:final_CHS_occupancies](data["tcrit"])
       end
   end
end

function TestIdealAFt(data)
    pyQmat = dcplikelihood.QMatrix(data["Q"] , data["nopen"])
    idealg = dcplikelihood.IdealG(pyQmat) 
    AFt = zeros(data["nopen"] , data["k"] - data["nopen"])
    x=ccall((:dcpIdealGXYt,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool ), AFt, data["Q"] , 2 , 5, data["t"][10], true )

    facts("Checking AF(t[10])") do
        @fact size(AFt) => (2 , 3)
        @fact AFt => idealg[:af](data["t"][10])
    end

end

function TestIdealFAt(data)
    pyQmat = dcplikelihood.QMatrix(data["Q"] , data["nopen"])
    idealg = dcplikelihood.IdealG(pyQmat) 

    FAt = zeros(data["k"] - data["nopen"], data["nopen"])
    x=ccall((:dcpIdealGXYt,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Csize_t,Csize_t,Float64,Bool ), FAt, data["Q"] , 2 , 5, data["t"][10], false )
    facts("Checking FA(t[10])") do
        @fact size(FAt) => (3 , 2)
        @fact FAt => idealg[:fa](data["t"][10])
    end
end
