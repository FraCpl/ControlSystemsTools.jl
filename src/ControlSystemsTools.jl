module ControlSystemsTools

using ControlSystems
using RobustAndOptimalControl
using GLMakie

export mag2db, db2mag, connectAuto
export nichols!, nicholsAxis

mag2db(x) = 20log10(x)
db2mag(x) = 10^(x/20)

function unwrapThis!(x, period=2π)      # Taken from JTools!
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k] - v,  y, RoundNearest)
	end
end

function nichols!(ax, sys, w; phaseBias::Int=0, kwargs...)
    mag, phase, w = bode(sys, w)
    lines!(ax, phase[:] .+ phaseBias*360.0, mag2db.(mag[:]); kwargs...)
end

function nichols!(ax, sys; phaseBias::Int=0, kwargs...)
    ub = ControlSystems.isdiscrete(sys) ? log10(π/sys.Ts) : 2
    w =  exp10.(range(-2; stop=ub, length=6000))
    nichols!(ax, sys, w; phaseBias=phaseBias, kwargs...)
end

function nicholsAxis(f; center::Int=-1)
    ct = 360*center + 180
    ax = Axis(f; limits=(ct - 180, ct + 180, -40, 40), xgridvisible=false, ygridvisible=false,
        xlabel="Open-Loop Phase [deg]", ylabel="Open-Loop Gain [dB]", title="Nichols Chart")
    phase = range(1, 359, 1000)*π/180;

    for mag in [6; 3; 2; 1; 0.5; 0; -0.5; -1; -2; -3; -6; -9; -12; -15; -20; -25; -30; -40; -50; -60]
        H = db2mag(mag).*exp.(im*phase)
        HH = H./(1 .- H)
        ph = atan.(imag.(HH), real.(HH))
        unwrapThis!(ph)
        mg = mag2db.(abs.(HH))

        lines!(ax, ph*180/pi .+ (ct - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
        lines!(ax, ph*180/pi .+ (ct + 360 - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
        lines!(ax, ph*180/pi .+ (ct - 360 - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
    end

    mag = -60:0.05:6
    for phs = [359; 355; 350:-10:10; 5; 1]
        H = db2mag.(mag).*exp(im*phs*pi/180)
        HH = H./(1 .- H)
        ph = atan.(imag.(HH), real.(HH))
        mg = mag2db.(abs.(HH))

        lines!(ax, ph*180/pi .+ 360.0.*(ph .≤ 0.0) .+ (ct - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
        lines!(ax, ph*180/pi .+ 360.0.*(ph .≤ 0.0) .+ (ct + 360 - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)
        lines!(ax, ph*180/pi .+ 360.0.*(ph .≤ 0.0) .+ (ct - 360 - 180), mg, linewidth=0.5, color=:grey, linestyle=:dash)

    end
    scatter!(ax, ct, 0, marker=:cross, color=:grey, markersize=8)
    return ax
end

connectAuto(SYS, in, out) = connectAuto(SYS)[out, in]
function connectAuto(SYS)
    # Extract all input and output symbols
    inSym = Symbol[]
    outSym = Symbol[]
    for S in SYS
        for s in S.u; push!(inSym, s); end
        for s in S.y; push!(outSym, s); end
    end
    if length(unique(outSym)) != length(outSym)
        @warn "Multiple systems are providing the same output."
    end

    # Build connections and identify system inputs
    connections = Pair{Symbol, Symbol}[]
    ins = Symbol[]
    for o in outSym
        if any(o .== inSym)
            push!(connections, o => o)
        end
    end
    for i in inSym
        if !any(i .== outSym)
            push!(ins, i)
        end
    end

    return connect(SYS, connections; w1=ins, unique=false)
end

end
