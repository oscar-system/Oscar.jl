using DelimitedFiles
using JSON
using Dates
using Statistics

injsonpath = "timing_summary.json"

indict = isfile(injsonpath) ? JSON.parsefile(injsonpath) : JSON.parse("{}")
filelist = readdir()
filelist = filter(endswith("csv"), filelist)
juliaVersion = join(split("$VERSION", ".")[1:2], ".")

if isempty(indict)
    indict["jobs"] = Dict()
end

for file in filelist
    (timestr, platform, juliaVersion, subset, commitHash) = split(file[1:end-4], "_")[2:end]
    # filenames on github must not contain colons so we convert to the correct timestamp format here
    timestamp = Dates.format(DateTime(timestr, dateformat"yyyy-mm-ddTHH-MM-SSZ"),
                             dateformat"yyyy-mm-ddTHH:MM:SS")
    commitAuthor = readchomp(`git show --no-patch --pretty=format:'%aN' $commitHash`)
    commitMessage = readchomp(`git show --no-patch --pretty=format:'%s' $commitHash`)
    if subset == "default"
        subset = "short"
    end
    testname = ":$platform: test $juliaVersion $subset"
    indata = readdlm(file, ',')

    totaltimesum = 0

    for entry in eachrow(indata)[2:end]
        ltestname = ":$platform: $juliaVersion $(entry[1])"
        timesum = sum(entry[2:end])
        totaltimesum = totaltimesum + timesum
        dict = Dict(
            [
                ("duration", timesum),
                ("date", timestamp),
                ("state", "passed"),
                ("commit", commitHash),
                ("message", commitMessage),
                ("author", commitAuthor)
            ]
        )
        if ! haskey(indict["jobs"], ltestname)
            indict["jobs"][ltestname] = Dict()
            indict["jobs"][ltestname]["recent"] = []
        end
        push!(indict["jobs"][ltestname]["recent"], dict)

        p = indict["jobs"][ltestname]["recent"]
        times = getindex.(p, "duration")
        if haskey(indict["jobs"][ltestname], "stats")
            t = indict["jobs"][ltestname]["stats"]
            t["count"] = t["count"] + 1
        else
            t = indict["jobs"][ltestname]["stats"] = Dict()
            t["count"] = 1
        end
        t["max_seconds"] = maximum(times)
        t["min_seconds"] = minimum(times)
        t["mean_seconds"] = mean(times)
        t["median_seconds"] = median(times)
        t["std_seconds"] = isnan(std(times)) ? 0 : std(times)
    end

    mydict = Dict(
        [
            ("duration", totaltimesum),
            ("date", timestamp),
            ("state", "passed"),
            ("commit", commitHash),
            ("message", commitMessage),
            ("author", commitAuthor)
        ]
    )
    if ! haskey(indict["jobs"], testname)
        indict["jobs"][testname] = Dict()
        indict["jobs"][testname]["recent"] = []
    end
    if haskey(indict["jobs"][testname], "stats")
        t = indict["jobs"][testname]["stats"]
        t["count"] = t["count"] + 1
    else
        t = indict["jobs"][testname]["stats"] = Dict()
        t["count"] = 1
    end
    push!(indict["jobs"][testname]["recent"], mydict)
    p = indict["jobs"][testname]["recent"]
    times = getindex.(p, "duration")
    t = indict["jobs"][testname]["stats"]
    t["count"] = t["count"] + 1
    t["max_seconds"] = maximum(times)
    t["min_seconds"] = minimum(times)
    t["mean_seconds"] = mean(times)
    t["median_seconds"] = median(times)
    t["std_seconds"] = isnan(std(times)) ? 0 : std(times)

    Base.Filesystem.rename(file, "$(file).done")
end

indict["generated_at"] = now()

JSON.json(injsonpath, indict)

