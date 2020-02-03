#making something to draw from a given (recorded)
#distribution of ISIs, when no stimulus is present
#and inserting spikes to "simulate" a stimulus.

#Author: Josephine Jefferson
#Start date: 03 Feb 2020
using MAT
using Plots
using Distributions
gr()

frequency=12.8 #frequency of simulated stimulus, in Hz
T=20 #Time in seconds
sensePhase=pi/4 #phase at which spike will be added
number=2 #number of spikes added each time you reach phase
every_what=1 #inverse frequency of adding spike - eg 2 would mean you add
#one spike every 2 cycles


#if want to determine # cycles instead of recording time
#cycles=8
#T= (1/frequency)*cycles

file = matopen("C:/Users/josep/OneDrive/Desktop/simulPhaseLock/data/CB94_19_144_05092013_spont_mgp.mat")
spontSpikes = read(file, "spt") #spike times in s
distrib= diff(spontSpikes[:])

Δt =.0001 #Sampling rate (s)
t = range(0, T, step = Δt) #Array of timestamps
SpikeTimes=zeros(length(t))
nSpikes=0

sptAsRadians=Float64[] #to store simulated spiketimes as radians.
meanVsize=0
meanVangle=0

function plotSpontISIs()
    display(histogram(distrib[:], nbins=25,
        title="Inter-Spike Interval Distribution of Spontaneous Afferent Dischacharge",
        xaxis="Inter-Spike Interval (s)", yaxis="Frequency", legend=false))
end

function genBaseSpikeTrain()
    nextSpikeTime=0
    while nextSpikeTime<T
        nextSpikeTime+=rand(distrib)
        global nSpikes+=1
        global SpikeTimes[nSpikes]=nextSpikeTime
    end
    global SpikeTimes[nSpikes]=0
    global nSpikes-=1
end

#insert spikes at a particular phase to simulate neuron being sensitive to
#that phase.
function insertSpikes(phase=pi/4, n=1, ev=1)
    #add in spikes at a particular phase
    #make sure to incrememnt nSpikes
    phlen=1/frequency #time/phase_length = radians/2pi
    firstInducedSpike = phase*(phlen/2pi)
    nextSpike = firstInducedSpike
    while nextSpike<T
        count=0
        while count<n
            global nSpikes+=1
            global SpikeTimes[nSpikes]=nextSpike
            count+=1
        end
        nextSpike+=phlen*ev
    end
end

function trimSpikeTimes()
    global SpikeTimes=SpikeTimes[1:nSpikes]
    sort!(SpikeTimes)
end

function plotSimISIs()
    display(histogram(diff(SpikeTimes[:]), nbins=25,
        title="Simulated Inter-Spike Interval Distribution of Spontaneous Afferent Dischacharge",
        xaxis="Inter-Spike Interval (s)", yaxis="Frequency", legend=false))
end

function convertSpt2Rad()
    phLen=1/frequency #time/phase_length = radians/2pi
    count=1
    sptAsRad=zeros(length(SpikeTimes))
    while count<=nSpikes
        answer=SpikeTimes[count]*(2pi/phLen)
        while answer>=2*pi
            answer=answer-(2*pi)
        end
        while answer<0
            answer+= 2*pi
        end
        sptAsRad[count]=answer
        count+=1
    end
    global sptAsRadians=sptAsRad[:]
end

#now we do stats!
function meanVectorCalc(print=false)
    n = length(sptAsRadians)
    sumCos=0.0
    for z in sptAsRadians[:]
        sumCos+=cos(z)
    end
    X=sumCos/n

    sumSin=0.0
    for z in sptAsRadians[:]
        sumSin+=sin(z)
    end
    Y=sumSin/n

    r=sqrt((X^2) + (Y^2))
    global meanVsize=r
    angle=atan(Y,X)
    global meanVangle = angle
    if print
        println("Vector strength: ", r)
        if angle<0
            angle+=2pi
        end
        println("Vector angle (rad): ", angle)
        println("Vector angle (deg): ", rad2deg(angle))
    end
end

#function displaying polarplot incl mean vector
function displayPolar()
    scale=1
    if meanVsize<0.01
        scale=0.02
    elseif meanVsize<0.05
        scale=0.1
    elseif meanVsize<0.1
        scale=0.2
    end
    x=zeros(length(sptAsRadians))
    count=1
    while count<=length(x)
        x[count]=scale+scale/10
        count+=1
    end

    polarplot=scatter(sptAsRadians[:], x, proj=:polar, label="Spikes")
    polarplot=plot!(0:meanVangle:meanVangle, 0:meanVsize:meanVsize, proj=:polar,
        yaxis=[0,scale+scale/5], linewidth=3, label="Mean Vector")
    display(polarplot)
end

function significanceTest(detail=false)
    spikes=sptAsRadians[:]
    r=meanVsize
    n=length(spikes)
    bigR=r*n
    rayZed=(bigR^2)/n
    if detail
        println("Rayleigh's R: ", bigR)
        println("Rayleigh's Z: ", rayZed)
    end
    P = exp(sqrt(1+4n+4(n^2-bigR^2))-(1+2n))
    println("p: ",P)
end



plotSpontISIs()
genBaseSpikeTrain()
insertSpikes(sensePhase, number, every_what)
trimSpikeTimes()
plotSimISIs()
convertSpt2Rad()
meanVectorCalc(true)
displayPolar()
significanceTest()
