using FactCheck
using PyCall
@pyimport dcprogs.likelihood as dcplikelihood

function TestIonGeometry()
    #setup code to match DCPROGS

    nopen = 2
    k = 5
    Q = [-3050        50  3000      0    0;
              2./3. -1502./3.     0    500    0;
              15       0 -2065     50 2000;
              0     15000  4000 -19000    0;
              0         0    10      0  -10;]

    burst_intervals = [0.1 0.2 0.1 0.2 0.15 0.16 0.18 0.05 0.1]
    burst_lengths = [3 1 5]
    tres = 5e-5
    tcrit = 1e-4
    useChs = true
    likelihood = [0.0]


    pybursts = cell(3)
    pybursts[1] = [0.1 , 0.2, 0.1]
    pybursts[2] = [0.2]
    pybursts[3] = [0.15, 0.16, 0.18, 0.05, 0.1]

    pylik = dcplikelihood.Log10Likelihood(pybursts, nopen, tres, tcrit)
    pyQmat = dcplikelihood.QMatrix(Q , nopen)
    
    x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int64,(Ptr{Float64},Ptr{Float64},Int64,Ptr{Int64},Csize_t,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,tcrit,useChs )
    pyll = pycall(pylik,PyAny,pyQmat)*log(10)


    facts("Check Likelihoods") do
        @fact likelihood[1] => pyll "dcprogs python gives $pyll"
        @fact x => 0
    end

    x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, 1e-5,tcrit,useChs )
    pylik = dcplikelihood.Log10Likelihood(pybursts, nopen, 1e-5, tcrit)
    pyll = pycall(pylik,PyAny,pyQmat)*log(10)

    facts("Testing tres change") do
        @fact likelihood[1] => pyll "dcprogs python gives $pyll"
        @fact x => 0
    end

    x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,1e-3,useChs )

    pylik = dcplikelihood.Log10Likelihood(pybursts, nopen, tres, 1e-3)
    pyll = pycall(pylik,PyAny,pyQmat)*log(10)

    facts("Testing tcrit change") do
        @fact likelihood[1] => pyll "dcprogs python gives $pyll"
        @fact x => 0
    end


    x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres, tcrit, 0 )


    pylik = dcplikelihood.Log10Likelihood(pybursts, nopen, tres, -tcrit)
    pyll = pycall(pylik,PyAny,pyQmat)*log(10)

    facts("CHS change") do
        @fact likelihood[1] => pyll "dcprogs python gives $pyll"
        @fact x => 0
    end


    burst_intervals = [0.1 0.2 0.1 0.2 0.15 0.16 0.18 0.05 0.1 0.10 0.18 0.02 0.1 0.03]
    burst_lengths = [3 1 5 5]

    x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,tcrit,useChs )

    push!(pybursts, [0.10, 0.18, 0.02, 0.1, 0.03])
    pylik = dcplikelihood.Log10Likelihood(pybursts, nopen, tres, tcrit)
    pyll = pycall(pylik,PyAny,pyQmat)*log(10)

    facts("Testing additional burst example from python dcprogs") do
        @fact likelihood[1] => pyll "dcprogs python gives $pyll"
        @fact x => 0
    end

    println("*** TIME TEST ***")
    println("\t1000 :ccall calls")
    @time for i=1:1000 x=ccall((:missed_events_likelihood,"libdcprogs_wrapper"),Int32,(Ptr{Cdouble},Ptr{Float64},Int64,Ptr{Int64},Int64,Ptr{Float64},Int64,Int64,Float64,Float64,Bool),likelihood,burst_intervals,length(burst_intervals),burst_lengths,length(burst_lengths),Q, nopen, k, tres,tcrit,useChs ) end

    println("\t1000 :PyCall calls")
    @time for i=1:1000 pylik = dcplikelihood.Log10Likelihood(pybursts, nopen, tres, tcrit); pyQmat = dcplikelihood.QMatrix(Q , nopen); pyll = pycall(pylik,PyAny,pyQmat)*log(10) end 


end
