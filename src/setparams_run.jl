using CSV
using DataFrames
using Glob

include("structs.jl")

directory_path = "/home/jakab/Desktop/Szakdoli/pisti_hiba_repr/1"
C_directory_path = "/home/jakab/Desktop/Szakdoli/pisti_hiba_repr/1/C"
text_files = glob("*.txt", directory_path)
C_files = glob("*.txt", C_directory_path)
data_collections = Dict()

function parse_content(content)
    lines = split(content, '\n', keepempty=false)
    # Function to safely parse a number from a string, removing commas
    safe_parse(s) = try parse(Float64, replace(s, "," => "")) catch; NaN; end
    
    # Check for a single column vector
    if all(line -> count(isspace, line) == 0, lines)
        return safe_parse.(lines)
        # Column vector
    elseif length(lines) == 1
        # Check for a single row vector
        return [safe_parse.(split(lines[1], isspace, keepempty=false))]
    else
        # It's a matrix
        return hcat([safe_parse.(split(line, isspace, keepempty=false)) for line in lines]...)'
    end
end

# Process each file
for file_path in text_files
    file_name_pre = replace(file_path, r"^/home/jakab/Desktop/Szakdoli/pisti_hiba_repr/1/self\." => "")
    file_name_pre2 = replace(file_name_pre, r"^/home/jakab/Desktop/Szakdoli/pisti_hiba_repr/1/" => "")
    file_name = replace(file_name_pre2, r"\.txt" => "")
    open(file_path, "r") do file
        content = read(file, String)
        data_collections[file_name] = parse_content(content)
    end
end

C = Array{Float64}(undef, 10, 10, 12)

for i in 1:12
    open(C_files[i], "r") do file
        content = read(file, String)
        parsed_array = parse_content(content)
        C[:,:,i] = parse_content(content)
    end
end


for key in keys(data_collections)
    println(key)
end

params = (12, 10, [], C, data_collections["D"], data_collections["W_ba"], vec(data_collections["species_n_reactions"]), vec(data_collections["species_n_splits"]), 
vec(data_collections["m"]), 0, 0, vec(data_collections["tau"]), vec(data_collections["alpha"]))

t_span = (0.0, 100)
prob = ODEProblem(equations2, vec(data_collections["y0"]), t_span, params)
solution = solve(prob, lsoda())

plot(solution)
savefig("test.png")






# Now data_collections contains all your data


