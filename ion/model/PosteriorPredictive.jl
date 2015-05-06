using PyCall
@pyimport dcprogs.likelihood as dcplikelihood

function AsymptoticRoots(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, tres::Float64)
   det = dcplikelihood.DeterminantEq(dcplikelihood.QMatrix(Q,nopen),tres)
    if !isOpen
        det = det[:transpose]()
    end

    #roots are returned as a tuple of the root and multiplicity. We just want the root 
    rootsandmulti = dcplikelihood.find_roots(det)
    
    roots = zeros(length(rootsandmulti))
    for i=1:length(roots)
        roots[i] = rootsandmulti[i][1]
    end
    return roots 
end

function IdealDistribution(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, t)

    idealg = dcplikelihood.IdealG(dcplikelihood.QMatrix(Q,nopen))

    if ndims(t) == 0
        t=[t] #make a vector
    end

    t = vec(t)
    density = zeros(length(t))
    if isOpen
        initial = idealg[:initial_occupancies] 
        G = idealg[:af](t) 
        final = ones(k - nopen)

        for i=1:length(t)
            density[i] = (initial * reshape(G[i,:,:],nopen,k-nopen) * final)[1]    
        end
    else
        initial = idealg[:final_occupancies]' #returned as a column vector so we need to transpose 
        G = idealg[:fa](t) 
        final = ones(nopen)
        for i=1:length(t)
            density[i] = (initial * reshape(G[i,:,:],k-nopen,nopen) * final)[1]    
        end
    end

    return t.*density 
    
end

function ExactDistribution(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, t, tres)
    if ndims(t) == 0
        t=[t]
    end

    exactg = dcplikelihood.MissedEventsG(dcplikelihood.QMatrix(Q,nopen), tres)

    t = vec(t)
    density = zeros(length(t))

    if isOpen
        initial = exactg[:initial_occupancies] 
        G = exactg[:af](t) 
        final = ones(k - nopen)

        for i=1:length(t)
            density[i] = (initial * reshape(G[i,:,:],nopen,k-nopen) * final)[1]    
        end
    else
        initial = exactg[:final_occupancies]' #returned as a column vector so we need to transpose 
        G = exactg[:fa](t) 
        final = ones(nopen)
        for i=1:length(t)
            density[i] = (initial * reshape(G[i,:,:],k-nopen,nopen) * final)[1]    
        end
    end

    return t.*density 

end

function ConditionalDistribution(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, tres::Float64, t, tlo::Float64, thi::Float64)
    if ndims(t) == 0
        t=[t]
    end

    dcpQ = dcplikelihood.QMatrix(Q,nopen)
    exactg = dcplikelihood.MissedEventsG(dcpQ, tres)
    conditionalDistribution = zeros(length(t))
    
    nclose = k-nopen
     
    KuHi = CumulativeExactSurvivor(Q, nopen, k, !isOpen, thi-tres , tres, false)
    KuLo = CumulativeExactSurvivor(Q, nopen, k, !isOpen, tlo-tres , tres, false)
    
    
    if isOpen
        phi = exactg[:final_occupancies]';
        cumulativeHi = (phi * KuHi * dcpQ[:fa] * expm(dcpQ[:aa]*tres) * ones(nopen))[1]
        cumulativeLo = (phi * KuLo * dcpQ[:fa] * expm(dcpQ[:aa]*tres) * ones(nopen))[1]
        for i=1:length(t)
            conditionalDistribution[i] = (phi * (KuHi - KuLo) * dcpQ[:fa] * expm(dcpQ[:aa]*tres) * exactg[:af](t[i]) * ones(nclose))[1]
            conditionalDistribution[i] = t[i]*conditionalDistribution[i]/(cumulativeHi-cumulativeLo)
        end    
    else
        phi = exactg[:initial_occupancies];
        cumulativeHi = (phi * KuHi * dcpQ[:af] * expm(dcpQ[:ff]*tres) * ones(nclose))[1]
        cumulativeLo = (phi * KuLo * dcpQ[:af] * expm(dcpQ[:ff]*tres) * ones(nclose))[1]
        for i=1:length(t)
            conditionalDistribution[i] = (phi * (KuHi - KuLo) * dcpQ[:af] * expm(dcpQ[:ff]*tres) * exactg[:fa](t[i]) * ones(nopen))[1]
            conditionalDistribution[i] = t[i]*conditionalDistribution[i]/(cumulativeHi-cumulativeLo)
        end    
    end
    return conditionalDistribution
end

function OpenConditionalMean(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, tres::Float64, tlo::Float64, thi::Float64)
    #E(t2 \mid [tlo < t1 < thi]) where t2 occurs after t1
    dcpQ = dcplikelihood.QMatrix(Q,nopen)
    exactg = dcplikelihood.MissedEventsG(dcpQ, tres)
    nclose = k - nopen
 
    CuHi = CumulativeExactDistribution(Q, nopen, k, !isOpen, thi , tres)[1] # returns an array
    CuLo = CumulativeExactDistribution(Q, nopen, k, !isOpen, tlo , tres)[1]

    KuHi = CumulativeExactSurvivor(Q, nopen, k , !isOpen, thi - tres, tres, false) #we want the cumulant K_{A|F}(u)
    KuLo = CumulativeExactSurvivor(Q, nopen, k , !isOpen, tlo - tres, tres, false)

    if isOpen
        conditionalMean = ((exactg[:final_occupancies]' * (KuHi - KuLo) * dcpQ[:fa]) * ((expm(dcpQ[:aa] * tres) *dARsds(dcpQ,nopen,nclose,tres) * dcpQ[:af]) *  expm(dcpQ[:ff] * tres) * ones(nclose)))[1]
    else
        error("Not implemented for shut sojourn preceeding open sojourn")
    end           

    return conditionalMean/(CuHi - CuLo)
end

function  dARsds (dcpQ::PyObject,nopen::Int64, nclose::Int64, tres::Float64)

    eFFt = expm(dcpQ[:ff]*tres);
    
    #P(sojourn in F states > tres)
    S_FF = eye(nclose,nclose) - eFFt #s=0

    #P(move to any F | start in any A) %s=0
    G_AF = (-dcpQ[:aa]^-1)*dcpQ[:af]

    #P(move to any state in A | start in any F) %s=0
    G_FA = (-dcpQ[:ff]^-1)*dcpQ[:fa]

    A1 = (tres * G_AF * (eFFt * G_FA)) - (G_AF * S_FF*(((dcpQ[:ff]^-1) * G_FA))) + ((-(dcpQ[:aa]^-1) * G_AF)*S_FF*G_FA)

    VA = eye(nopen,nopen) - ((G_AF*S_FF)*G_FA)
    A2 = dcpQ[:aa]^(-1)-(A1*(VA^-1));
    deriv = ((VA^-1) * A2) * dcpQ[:aa]^-1;

end

function UnivariateConditionalMean(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, tres::Float64, tlo::Float64, thi::Float64)
    #E(t \mid tlo < t < thi)
    dcpQ = dcplikelihood.QMatrix(Q,nopen)
    exactg = dcplikelihood.MissedEventsG(dcpQ, tres)
    nclose = k - nopen

    CuHi = CumulativeExactDistribution(Q, nopen, k, isOpen, thi , tres)
    CuLo = CumulativeExactDistribution(Q, nopen, k, isOpen, tlo , tres)

    KvHi = CumulativeExactSurvivor(Q, nopen, k , isOpen, thi - tres, tres, true) #we want the mean cumulant vK_{A}(t)
    KvLo = CumulativeExactSurvivor(Q, nopen, k , isOpen, tlo - tres, tres, true)
    
    
    if isOpen
        phi = exactg[:initial_occupancies]
        conditionalMean =  (phi * (KvHi - KvLo) *dcpQ[:af] * (expm(dcpQ[:ff]*tres) * ones(nclose))/(CuHi-CuLo))[1]
    else
        phi = exactg[:final_occupancies]
        conditionalMean =  (phi' * (KvHi - KvLo) * dcpQ[:fa] * (expm(dcpQ[:aa]*tres) * ones(nopen))/(CuHi-CuLo))[1]

    end

    return tres + conditionalMean;

end

function CumulativeExactDistribution(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, t, tres::Float64)
    if ndims(t) == 0
        t=[t]
    end

    dcpQ = dcplikelihood.QMatrix(Q,nopen)
    exactg = dcplikelihood.MissedEventsG(dcpQ, tres)

    t = vec(t)
    density = zeros(length(t))

    if isOpen
        initial = exactg[:initial_occupancies] 
        final = ones(k - nopen)

        for i=1:length(t)
            density[i] = (initial * CumulativeExactSurvivor(Q , nopen, k, isOpen, t[i] - tres, tres, false) * dcpQ[:af]  * expm(dcpQ[:ff]*tres) * final)[1]    
        end
    else
        initial = exactg[:final_occupancies]' #returned as a column vector so we need to transpose 
        final = ones(nopen)
        for i=1:length(t)
            density[i] = (initial * CumulativeExactSurvivor(Q , nopen, k, isOpen, t[i]-tres, tres, false) * dcpQ[:fa] * expm(dcpQ[:aa] * tres) * final)[1]    
        end
    end
    return density
end

function CumulativeExactSurvivor(Q::Array{Float64,2}, nopen::Int64, k::Int64, isOpen::Bool, u::Float64, tres::Float64, meanCumulSurvivor::Bool)

    if isOpen
        Ku = zeros(nopen,nopen)
    else
        Ku = zeros(k - nopen, k - nopen)
    end

    lambdas = eigvals(-Q)

    if u <= 0
        return Ku
    end

    survivor = dcplikelihood.ExactSurvivor(dcplikelihood.QMatrix(Q,nopen),tres)

    #we want either the exact cumulant upto u if u < 2 * tres or the whole cumulant to be added to by the approx solution
    if u > 2 * tres
        v= 2 * tres
    else
        v = u
    end

    for i=0:length(lambdas)-1
        for m = 0:convert(Int64,floor(v/tres))
            for r = 0:m
                if meanCumulSurvivor
                    cumulant = ((tres* m * iCG(v-(m*tres),lambdas[i+1],r)) + iCG(v-(m*tres),lambdas[i+1],r+1))
                else
                    cumulant = iCG(v-(m*tres),lambdas[i+1],r)
                end

                if isOpen
                    Ku+=((-1.0)^m) * survivor[:recursion_af](i,m,r) * cumulant#iCG(v-(m*tres),lambdas[i+1],r)
                else
                    Ku+=((-1.0)^m) * survivor[:recursion_fa](i,m,r) * cumulant#iCG(v-(m*tres),lambdas[i+1],r)
                end
            end
        end
    end

    if u > (2 * tres)
        #need to add additional cumulation from approximate survivor function
        approxsurvivor = dcplikelihood.ApproxSurvivor(dcplikelihood.QMatrix(Q,nopen),tres)
        if isOpen
            components = approxsurvivor[:af_components]
        else
            components = approxsurvivor[:fa_components]
        end
        if meanCumulSurvivor
            for i=1:length(components)
                timeConstant = -1/components[i][2]
                Ku+=(components[i][1] * timeConstant * (((timeConstant+ (2*tres)) * (exp((-2*tres)/timeConstant))) - ((timeConstant+u) * exp(-u/timeConstant))))
            end
        else
            for i=1:length(components)
               timeConstant = -1/components[i][2];
               Ku=Ku+(components[i][1]*timeConstant*(exp((-2*tres)/timeConstant) - exp(-u/timeConstant)));            end
        end
    end
    return Ku
end

function iCG(u::Float64, lambda::Float64, r::Int64)
    #Evaluates the integral \int^u_0 v^rexp(-lambda*v) dv
    if lambda < 1e-16
        icg = u^(r+1)/(r+1);
    else
        summate = 0;
        for j=0:r
            summate = summate + exp(-lambda*u) * ((lambda*u)^j)/factorial(j);
        end
        summate = 1 - summate;
        icg = summate * (factorial(r)/lambda^(r+1));
    end
    return icg
end

