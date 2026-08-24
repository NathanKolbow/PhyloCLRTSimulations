using Dates

function format_elapsed(seconds::Real)
    # Convert total seconds to Millisecond for Dates rounding
    ms = Millisecond(round(Int, seconds * 1000))
    
    # Very fast calls (< 1 second): show decimal seconds or ms directly
    if ms < Millisecond(1000)
        return "$(round(seconds * 1000, digits=2)) ms"
    end
    
    # Use Dates.CanonicalDecode to break down ms into Period objects (e.g., [1 Minute, 12 Seconds])
    period = Dates.canonicalize(ms)
    
    # Format periods neatly (e.g., "1m 12s", "2h 15m 3s")
    parts = String[]
    for p in period.periods
        if p isa Hour
            push!(parts, "$(p.value)h")
        elseif p isa Minute
            push!(parts, "$(p.value)m")
        elseif p isa Second
            push!(parts, "$(p.value)s")
        elseif p isa Millisecond
            push!(parts, "$(p.value)ms")
        end
    end
    
    return join(parts, " ")
end

macro simlog(stringarg, fxncall)
    return quote
        print($(esc(stringarg)))
        flush(stdout)
        
        t0 = time_ns()
        val = $(esc(fxncall))
        elapsed = (time_ns() - t0) / 1e9
        
        # Call the helper function to format the output
        println(" (elapsed: ", format_elapsed(elapsed), ")")
        
        val
    end
end

macro simlog(stringarg)
    return quote
        println($(esc(stringarg)))
        flush(stdout)
    end
end

macro simlog_assertion(stringarg, assertion)
    return quote
        print($(esc(stringarg)))
        flush(stdout)

        if $(esc(assertion))
            printstyled(" TRUE\n", color=:green)
        else
            printstyled(" FALSE\n", color=:red)
        end
    end
end
