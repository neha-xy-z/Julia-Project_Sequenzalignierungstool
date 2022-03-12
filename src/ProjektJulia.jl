module ProjektJulia

# Needed Modules
using FASTX
using ZipFile
import FASTX.FASTQ
import Base

# globals
global adapter_seq = Ref("")
global adapter_seq_str = Ref("")
global fastqFile = Ref("")

#julia> pwd()
#julia> cd("C:\\Users\\PC\\ProjektJulia\\src")
#julia> include("ProjektJulia.jl")
# C:\Users\PC\Downloads\SRR14608302\SRR14608302.fastq

# Needed file
function askForFileAndDoesFileExist()
    print("\nEnter path of your FASTQ file. \n\n")
    f = readline()
    if isfile(f)
        if endswith(f, ".fastq")
            return f
        else
            print("File input is not a FASTQ File. Check your path again.")
        end
    else
        print("File not found. Check your path again and enter on request. Or type 'q' to quit.\n\n")
        if readline() == "q"
            print("Bye.")
            return
        else 
            askForFileAndDoesFileExist()
        end
    end
end

fastqFile = askForFileAndDoesFileExist()

# Set Path
function set_path(file_path::String)
    if Sys.islinux()
        file_fastx = string(split(file_path, "/")[end])
        file_name = string(split(file_fastx, ".")[1])
        work_path = string(replace(file_path, "/"*file_fastx => ""))
    else
        file_fastx = string(split(file_path, "\\")[end])
        file_name = string(split(file_fastx, ".")[1])
        work_path = string(replace(file_path, "\\"*file_fastx => ""))
    end
    return file_fastx, file_name, work_path
end

file_fastx, file_name, work_path = set_path(fastqFile)
    
# Set Path to used Softwares    
function setPathToSoftware()
    curr_dir = pwd()
    if Sys.islinux()
        curr_dir = replace(curr_dir, "/src" => "")

        path_to_cutadapt = curr_dir*"/cutadapt"
        path_to_fastqc = curr_dir*"/fastqc/FastQC"
        path_to_vsearch = curr_dir*"/vsearch/vsearch/bin"
    else
        curr_dir = replace(curr_dir, "\\src" => "")

        path_to_cutadapt = curr_dir*"\\cutadapt"
        path_to_fastqc = curr_dir*"\\fastqc\\FastQC"
        path_to_vsearch = curr_dir*"\\vsearch\\vsearch\\bin"
    end
    return path_to_cutadapt, path_to_fastqc, path_to_vsearch
end

path_to_cutadapt, path_to_fastqc, path_to_vsearch = setPathToSoftware()

# Cutadapt Run for different OS Linux, Windows
function cutadapt(adapter_seqs_str::String, file_fastx::String)
    if Sys.islinux()
        run(`cutadapt -j 0 $adapter_seqs_str -o without_adapters.fastq $file_fastx`)
    else
        cd(path_to_cutadapt)
        run(`cutadapt-3.4.exe $adapter_seqs_str -o without_adapters.fastq $file_fastx`) # Windows
        curr_dir = replace(path_to_cutadapt, "\\cutadapt" => "")
        cd(curr_dir*"\\src")
    end
end


# FastQC and Cutadapt Run for different OS Linux, Windows
function func()
    print("\n\nDo you want to run FastQC to get a possible Source for an Adapter Sequence, which could be found in your input-file?\n If yes, type 'y'. If no, type 'n'. You can enter an individual Adapter sequence afterwards on request.\n\n")
    
    answer = readline()
    if answer == "y"
        if Sys.islinux()
            # FastQC Run for OS Linux

            function fastqcUnix(input_path::String, output_path::String)
                run(`fastqc $input_path -o $output_path`)
            end
            fastqcUnix(fastqFile, path_to_fastqc)
            
        else 
            # FastQC Run for Windows
            
            function fastqcWindows(input_path::String, output_path::String)
                cd(path_to_fastqc)
                run(`java -Xmx250m -classpath .";"./sam-1.103.jar";"./jbzip2-0.9.jar uk.ac.babraham.FastQC.FastQCApplication $input_path`)
                curr_dir = replace(path_to_fastqc, "\\fastqc\\FastQC" => "")
                cd(curr_dir*"\\src")
            end
            
            #fastqcWindows(fastqFile, path_to_fastqc)
        end

        print("\nFastQC successful.\n")
        print("\nReading adapter sequences...\n")

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

            for adapter_seq in seqs
                adapter_seq_str *= "-a $adapter_seq "
            end

            return seqs, sources, string(strip(adapter_seq_str))
        end

        fastqc_zip_file = replace(fastqFile, ".fastq" => "_fastqc.zip")
        print(fastqc_zip_file)
        adapter_seqs, adapter_sources, adapter_seq_str = read_adapter(fastqc_zip_file)
        print("\nSuccessful.\n")

        println("Adapters found:")
        println(adapter_seqs, adapter_sources, adapter_seq_str)

        # Run Cutadapt
        cutadapt(adapter_seq_str, fastqFile)
   
    elseif answer=="n"
        print("\nPlease type in your adapter sequence:\n")
        adapter_seq = readline()
        cutadapt(adapter_seq, fastqFile)
    else
        print("Wrong input.")
        func()
    end
end

# call function func for FastQC and Cutadapt Run
# TODO change name
func()

# Read adapterfree fastq File and store into records
function readFastq(file_path::String)
    print("Reading file without adapters...")
    reader = FASTQ.Reader(open(file_path, "r"))
    records = []
    for record in reader
        push!(records, record)
    end
    close(reader)
    return records
end

# find file without adapters 
file_without_adapters_file = path_to_cutadapt*"\\file_without_adapters.fastq"

# read fastq into records
fastq_records = readFastq(file_without_adapters_file)

# Read adapterfree fastq File and store into records
function readFasta(file_path::String)
    reader = FASTA.Reader(open(file_path, "r"))
    records = []
    for record in reader
        push!(records, record)
    end
    close(reader)
    return records
end

# Collect garbage and free memory 
function free_memory(variable)
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
threshold_Param = getParamForTrimming()

# Trimming
function trimming(fileName::String, org_records::Vector{Any}, threshold::UInt8, min_length::Int)
    
    println("\n Preparing trimming...")
    fileName = replace(fileName, ".fastq" => "")
    fn = string(fileName)*"_trimmed.fastq"
    file = FASTQ.Writer(open(fn, "w"),true)
    
    for record in org_records
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

    # return trimmed_records
    close(file)
    free_memory(fastq_records)
    println(string(fn))
    println("Successful.")
    return string(fn)
end

#TODO threshold and min_length should als parameter 
fileForVSearch = trimming(fastqFile, fastq_records, threshold_Param, 30) 

# VSearch Run for different OS Linux, Windows
if Sys.islinux()
    function vsearchrun(input_path::String)
        run(`vsearch --uchime3_denovo $input_path --nonchimeras nonchimeras.fasta`)
    end
else
    println(pwd())
    
    function vsearchrun(input_path::String)
        cd(path_to_vsearch)
        run(`vsearch --uchime3_denovo $input_path --nonchimeras nonchimeras.fasta`)
        curr_dir = replace(path_to_vsearch, "\\vsearch\\vsearch\\bin" => "")
        cd(curr_dir*"\\src")
    end
end

# run vsearch
vsearchrun(fileForVSearch)

# find vsearch file and read into records
vsearch_file = path_to_vsearch*"\\nonchimeras.fasta"
recs = readFasta(vsearch_file)

# Enter kmer size
function foo()
    
    println("\nEnter a kmer-size:\n")
    n = parse(UInt8, readline())
    if n != UInt8
        println("\nWrong input.\n")
        foo()
    else
        return n
    end

end

# Save size into variable for func newReads()
kmer_size = foo()


# Count frequent kmers
function kmers_Histo(k, records)
    kmer_dict = Dict{String,Integer}()

    for record in records
        rec_seq = string(FASTQ.sequence(record))

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

        seq = string(FASTQ.sequence(rec))
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

newReads(kmers_Histo_Dict, recs, threshold_Param, kmer_size)

end # module