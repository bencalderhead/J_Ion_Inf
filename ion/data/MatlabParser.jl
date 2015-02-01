function MatlabParser(datafile)
    file = matopen(datafile)
    #TODO need to handle invalid files here

    bursts = read(file,"bursts")
    concs = read(file,"concs")
    tcrit = read(file,"tcrit")
    tres = read(file,"tres")
    useChs =  read(file,"useChs")
    close(file)
 
    experimentno = length(concs)
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
    end
    pybursts = cell(experimentno)
    for i=1:experimentno
        set = cell(length(bursts[i]))
        for j=1:length(bursts[i])
            if length(bursts[i][j]) > 1
                set[j] = vec(bursts[i][j])
            else
                set[j] = [bursts[i][j]]
            end
        end
        pybursts[i] = set
    end

    return {"pybursts" => pybursts, "experiment_nos" => experimentno, "burst_intervals" => bsts , "burst_lengths" => burst_lengths, "burst_nos" => burstnos, "concs"=> concs, "tcrit" => tcrit, "tres" => tres, "useChs" => useChs}
end        
            
