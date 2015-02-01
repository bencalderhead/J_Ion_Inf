using PyPlot
using IonTest

include("/home/michaelepstein/J_Ion_Inf_Buggy/ion/model/PosteriorPredictive.jl")

#open distributions
isOpen=true
openr = AsymptoticRoots(IonTest.Q, IonTest.nopen ,IonTest.k ,isOpen,IonTest.tres)
t=logspace(log10(0.00001),log10(20*(-1/maximum(openr))),512)

idealOpen=IdealDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, t)

semilogx(t, idealOpen, "r--")
for i in [2.5e-5,5e-5,1e-4,2e-4] 
    exactOpen = ExactDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, t, i)
    semilogx(t,exactOpen, "b--")
end
title("Unconditional Open Times")

shutr = AsymptoticRoots(IonTest.Q, IonTest.nopen ,IonTest.k ,!isOpen,IonTest.tres)
t=logspace(log10(0.00001),log10(20*(-1/maximum(shutr))),512)
idealShut = IdealDistribution(IonTest.Q, IonTest.nopen, IonTest.k, !isOpen, t)

figure()
plot(t, idealShut, "r--")
for i in [2.5e-5,5e-5,1e-4,2e-4] 
    exactShut = ExactDistribution(IonTest.Q, IonTest.nopen, IonTest.k, !isOpen, t, i)
    semilogx(t,exactShut, "b--")
end
title("Unconditional Shut Times")

#conditional distributions
t=logspace(log10(0.00001),log10(20*(-1/maximum(openr))),512)
figure()
title("Conditional subplots Figure 8a and b from 1996")
exactOpen = ExactDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, t, 5e-5)
subplot(1,2,1)
plot(t, exactOpen, "r--")
conditionalOpen = ConditionalDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, 5e-5, t, 5e-5,1.5e-4)
semilogx(t,conditionalOpen, "b--")

subplot(1,2,2)
plot(t,exactOpen,"r--")
conditionalOpen = ConditionalDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, 5e-5, t, 1e-2,1.0)
semilogx(t,conditionalOpen, "b--")


figure()
title("Conditional subplots Figure 8c and d from 1996")
exactOpen = ExactDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, t, 2e-4)
subplot(1,2,1)
plot(t, exactOpen, "r--")
conditionalOpen = ConditionalDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, 2e-4, t, 2e-4,1e-3)
semilogx(t,conditionalOpen, "b--")

subplot(1,2,2)
plot(t,exactOpen,"r--")
conditionalOpen = ConditionalDistribution(IonTest.Q, IonTest.nopen, IonTest.k, isOpen, 2e-4, t, 1e-2,1.0)
semilogx(t,conditionalOpen, "b--")

#mean continuous conditional distribution

ranges = [(5e-5,1e-4),(1e-4,2e-4),(2e-4,1e-3),(1e-3,3e-1),(3e-1,3.0)]
conditionalOpenMean = zeros(length(ranges))
preceedingMeanShut = zeros(length(ranges))

for i=1:4
    conditionalOpenMean[i] = OpenConditionalMean(Q, 2, 5, true, 5e-5 , ranges[i][1],ranges[i][2])
    preceedingMeanShut[i] = UnivariateConditionalMean(Q, 2, 5, false, 5e-5 , ranges[i][1],ranges[i][2])
end
figure()
scatter(preceedingMeanShut,conditionalOpenMean)




