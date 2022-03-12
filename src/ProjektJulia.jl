# Needed Modules
using FASTX
using ZipFile

# Set Path For Linux
function get_path_variable(file_path::String)

    if Sys.islinux()
        file_fastx = string(split(file_path, "/")[end])     # xxx.fastq
        file_name = string(split(file_fastx, ".")[1])       # xxx
        work_path = string(replace(file_path, "/"*file_fastx => ""))    # xx/xx/xx
    else
        file_fastx = string(split(file_path, "\\")[end])
        file_name = string(split(file_fastx, ".")[1])
        work_path = string(replace(file_path, "\\"*file_fastx => ""))
    end
    return file_name, work_path
end

# Set Path to used Softwares for Windows   
function setPathToSoftwareForWindows()
    curr_dir = pwd()
    curr_dir = replace(curr_dir, "\\src" => "")

    path_to_cutadapt = curr_dir*"\\cutadapt"
    path_to_fastqc = curr_dir*"\\fastqc\\FastQC"
    path_to_vsearch = curr_dir*"\\vsearch\\vsearch\\bin"
   
    return path_to_cutadapt, path_to_fastqc, path_to_vsearch
end

path_to_cutadapt, path_to_fastqc, path_to_vsearch = setPathToSoftwareForWindows()

# FastQC Unix Function
function fastqcUnix(input_file::String)
    run(`fastqc $input_file`)
end

# FastQC Memory Function
function fastqcWindows(input_fastqc_path::String, input_path::String, output_path::String)
    cd(input_fastqc_path)
    run(`java -Xmx250m -classpath .";"./sam-1.103.jar";"./jbzip2-0.9.jar uk.ac.babraham.FastQC.FastQCApplication $input_path`)
    curr_dir = replace(input_fastqc_path, "\\fastqc\\FastQC" => "")
    cd(curr_dir*"\\src")
end

# Cutadapt Run for Linux
function cutadaptUnix(adapter_seqs_str::String, input_file::String)
    run(`cutadapt -j 0 $adapter_seqs_str -o without_adapters.fastq $input_file`)
end

# Cutadapt Run for Windows
function cutadaptWindows(input_cutadapt_path::String, adapter_seqs_str::String, file_fastx::String)
    cd(input_cutadapt_path)
    run(`cutadapt-3.4.exe $adapter_seqs_str -o without_adapters.fastq $file_fastx`) # Windows
    curr_dir = replace(input_cutadapt_path, "\\cutadapt" => "")
    cd(curr_dir*"\\src")
end 

# Function for reading adapter from FastQC File
function read_adapter(fastqc_result_path::String)
            
    global fastqc_res = ZipFile.Reader(fastqc_result_path)

    flag = false
    seqs = String[]
    counts = String[]
    percentages = String[]
    sources = String[]
    adapter_seq_str = ""

    for file in fastqc_res.files
        if contains(file.name, "fastqc_data.txt")
            for line in eachline(file)
                if contains(line, ">>Overrepresented")
                    flag = true
                    continue
                end

                if flag == true && contains(line, ">>END_MODULE")
                    flag = false
                    break
                end

                if flag && !contains(line, "#")
                    seq, count, percentage, source = split(line, "\t")
                    if string(source) != "No Hit"
                        push!(seqs, seq)
                        push!(counts, count)
                        push!(percentages, percentage)
                        push!(sources, source)
                    end
                end
            end
        end
    end 

    if length(seqs) == 0
        #TODO
        adapter_seq_str = "no adapter seq"
    else

        for adapter_seq in seqs
            adapter_seq_str *= "-a $adapter_seq "
        end
    end

    return seqs, sources, string(strip(adapter_seq_str))
end

# FastQC and Cutadapt Run for different OS Linux, Windows
function func()
    org_fastq_file_path = ""
    work_path = ""
    file_name = ""

    while true
        print("\nEnter path of your FASTQ file. Or type 'q' to quit.\n>>> ")

        input = string(strip(readline()))
        if input == "q"
            println("Bye.")
            exit()
        elseif isfile(input)
            if endswith(input, ".fastq")
                org_fastq_file_path = input
                break
            else
                println("ERROR:\nFile input is not a FASTQ File. Check your path again.")
            end
        else
            println("ERROR:\nFile not found. Check your path again and enter on request.")
        end
    end

    # 
    file_name, work_path = get_path_variable(org_fastq_file_path)

    cd(work_path)

    while true 
        print("\nDo you want to run FastQC to get a possible Source for an Adapter Sequence, which could be found in your input-file?
                 \nIf yes, type 'y'.
                 \nIf no, type 'n'. Then you can enter an individual Adapter sequence afterwards on request.
                 \nOr type 'q' to quit.
                 \n>>> ")

        input = readline()

        if input == "q"
            println("Bye.")
            exit()

        elseif input == "y"
        
            if Sys.islinux()
                # FastQC Run for OS Linux
                fastqcUnix(file_name*".fastq")
            else 
                # FastQC Run for Windows
                fastqcWindows(path_to_fastqc, org_fastq_file_path, path_to_fastqc)
            end

            println("\nFastQC successful.")

            println("\nReading adapter sequences...")

            if Sys.islinux()
                adapter_seqs, adapter_sources, adapter_seq_str = read_adapter(file_name*"_fastqc.zip")
            else
                adapter_seqs, adapter_sources, adapter_seq_str = read_adapter(work_path*"\\"*file_name*"_fastqc.zip")
            end
            
            println("Adapters found:")

            if startswith(adapter_seq_str, "no")
                println("\t"*adapter_seq_str)
                break
            else
                for i = 1 : length(adapter_seqs)
                    println(adapter_seqs[i], adapter_sources[i])
                end

                # Run Cutadapt
                if Sys.islinux()
                    cutadaptUnix(adapter_seq_str, file_name*".fastq")
                    break
                else
                    #  Windows
                    cutadaptWindows(path_to_cutadapt, adapter_seq_str, work_path*"\\"*file_name*".fastq")
                    break
                end
            end

        elseif input == "n"
            print("\nPlease type in your adapter sequence:\n>>> ")

            adapter_seq = readline()

            if contains(adapter_seq, r"[^atcgATCG]")
                println("\n warning check again")
            else
                # Run Cutadapt
                if Sys.islinux()
                    cutadaptUnix("-a "*adapter_seq, file_name*".fastq")
                    break
                else
                    cutadaptWindows(path_to_cutadapt, adapter_seq_str, work_path*"\\"*file_name*".fastq")
                    break
                end
            end

        else
            println("ERROR:\nWrong input.")
        end
    end
end
func()

# Read adapterfree fastq File and store into records
function readFastq(input_file::String)
    # println("Reading file without adapters...")
    reader = FASTQ.Reader(open(input_file, "r"))
    records = []
    for record in reader
        push!(records, record)
    end
    close(reader)
    return records
end

# Read adapterfree fastq File and store into records
function readFasta(input_file::String)
    reader = FASTA.Reader(open(input_file, "r"))
    records = []
    for record in reader
        push!(records, record)
    end
    close(reader)
    return records
end

# Collect garbage and free memory 
function free_memory(variable::Any)
    empty!(variable)
    GC.gc()
end

# Get threshold parameter for the trimming function
function getParamForTrimming()
    print("\nEnter threshold for trimming (3/4/5): \n\n")
    n = parse(UInt8, readline())
    if n==3 || n==4 || n==5
        return n
    else
        print("You can only enter a threshold size of 3, 4 or 5. Please enter again.")
        getParamForTrimming()
    end
end

# find file without adapters 
file_without_adapters_file = path_to_cutadapt*"\\without_adapters.fastq"

# read fastq into records
fastq_records = readFastq(file_without_adapters_file)

threshold_Param = getParamForTrimming()

# Trimming
function trimming(input_file::String, input_records::Vector{Any}, threshold::UInt8, min_read_length::Int)
    
    println("\n Preparing trimming...\n")
    fileName = replace(fileName, ".fastq" => "")
    fn = string(fileName)*"_trimmed.fastq"
    file = FASTQ.Writer(open("trimmed.fastq", "w"), true)
    
    for record in input_records
        quality_list = quality(record)
        seq_list = sequence(String, record)
        count = 0
        position = 0

        # 'N' und short-reads skip
        if length(seq_list) < min_length || contains(seq_list, 'N')
            continue
        end
        
        for base_index = 1 : length(quality_list)    
            if Int64(quality_list[base_index]) <= 20
                count += 1
                if count == threshold
                    position = base_index
                    break
                end
            end
        end

        # low quality reads 
        if count == threshold 

            # only work with long-reads 
            if (position-threshold) >= min_length

                # New Sequence with new Quality
                new_seq = sequence(String, record, (1:position-threshold))
                new_quality = quality(record, 33, (1:position-threshold))
                # Get the same identifiers and description from the old Record and give it to the new Record
                ident = identifier(record)
                desc = description(record)
                new_record = FASTQ.Record(ident, desc, new_seq, new_quality)
                write(file, new_record)
            end
        else
            write(file, record) 
        end
    end
    return trimmed_records
    close(file)
    free_memory(fastq_records)
    # println(string(fn))
    println("Successful.")
    return string(fn)
end

#TODO threshold and min_length should als parameter 
fileForVSearch = trimming(work_path*"\\"*file_name*".fastq", fastq_records, threshold_Param, 30) 

# VSearch Run for different OS Linux, Windows
function vsearchrun(input_vsearch_path::String, input_path::String)
    if Sys.islinux()
        run(`vsearch --uchime3_denovo $input_path --nonchimeras nonchimeras.fasta`)
    else
        # println(pwd())
        cd(input_fastqc_path)
        run(`vsearch --uchime3_denovo $input_path --nonchimeras nonchimeras.fasta`)
        curr_dir = replace(input_fastqc_path, "\\vsearch\\vsearch\\bin" => "")
        cd(curr_dir*"\\src")
    end
end

# run vsearch
vsearchrun(path_to_vsearch, fileForVSearch)

# # find vsearch file and read into records
vsearch_file = path_to_vsearch*"\\nonchimeras.fasta"
recs = readFasta(vsearch_file)

function kmer_input()
    
    println("\nEnter a kmer-size:\n")
    n = parse(UInt8, readline())
    if n != UInt8
        println("\nWrong input.\n")
        kmer_input()
    else
        return n
    end
end

kmer_size = kmer_input()

# Write FASTQ File filled with record
# function writeFASTQ(fileName::String, completely_Filtered::Vector{Any})
#     file = FASTQ.Writer(open(string(fileName)*".fastq", "w"), true)
#     for record in completely_Filtered
#         write(file, record)
#     end
#     close(file)
# end

# Count frequent kmers
function kmers_Histo(k, records)
    kmer_dict = Dict{String,Integer}()

    for record in records
        rec_seq = string(FASTA.sequence(record))

        for i = 1:((length(rec_seq) - k) + 1)
            sequence_k = rec_seq[i:((i+k)-1)]
            if !haskey(kmer_dict, sequence_k)
                kmer_dict[sequence_k] = 1
            else
                kmer_dict[sequence_k] += 1
            end
        end
    end
    return kmer_dict
end

kmers_Histo_Dict = kmers_Histo(kmer_size, recs)

# Find neighbours from a kmer
function hammingDistance(kmer)

    alphabet = ['A','C','G','T']
    neighbors = []

    for all_pos = reverse(1:length(kmer))

        for character in alphabet

            new_kmer = collect(kmer)

            if character != new_kmer[all_pos]
                new_kmer[all_pos] = character
                push!(neighbors, join(new_kmer))
            end
        end

    end

    return neighbors
end

# Correct Reads into free-error Reads
function newReads(kmersHisto, recs, t, k)

    file = FASTA.Writer(open("new_errorfree_reads.fasta", "w"), 1000)
    
    for rec in recs

        seq = string(FASTA.sequence(rec))
        n_Kmers = (length(seq) - k) + 1

        for i = 1:n_Kmers

            kmer = seq[i:(i+k)-1]

            if haskey(kmersHisto, kmer)

                if kmersHisto[kmer] < t

                    all_neighbors = hammingDistance(kmer)

                    for nk in all_neighbors
                        if haskey(kmersHisto, nk)

                            if kmersHisto[nk] >= t

                                first_part = seq[1:i-1]
                                kmer_part = nk
                                second_part = seq[k+i:end]
    
                                new_seq = string(first_part, kmer_part, second_part)
                                id = identifier(rec)
    
                                new_record = FASTA.Record(id,new_seq)
                                write(file, new_record)

                                break
                            end
                        end
                    end
                    free_memory(all_neighbors)
                end
            end
        end
    end
    close(file)
end 

# newReads(kmers_Histo_Dict, recs, threshold_Param, kmer_size)

