module BTK
using LinearAlgebra
using Plots
using LaTeXStrings
using QuadGK
using MAT

function get_params(E)
    ϵk = √Complex(E^2 - Δ^2)
    uk = √(1/2*(1 + ϵk/E))
    vk = √(1/2*(1 - ϵk/E))
    # multiply num and den with c^2
    kplus = √(2*mc2*(μ + ϵk)/(ħ*c)^2)
    kminus = √(2*mc2*(μ - ϵk)/(ħ*c)^2)
    qplus = √(2*mc2*(μ + E)/(ħ*c)^2)
    qminus = √(2*mc2*(μ - E)/(ħ*c)^2)
    return (uk, vk, kplus, kminus, qplus, qminus)
end

function calc_amplitudes(E, H, params)
    uk, vk, kplus, kminus, qplus, qminus = params
    imthing = 2*mc2*H*1im/(ħ*c)^2
    coeff1 = imthing + qplus  
    coeff2 = imthing - qminus
    lhs = [0 1 -uk -vk;
           1 0 -vk -uk;
           0 coeff1 kplus*uk -kminus*vk;
           coeff2 0 kplus*vk -kminus*uk]
    rhs = [-1; 0; qplus-imthing; 0]
    sol = lhs\rhs
    return sol 
end

function vel(k)
    return ħ*c^2*k/mc2
end

function coeffs(amps, params, all)
    uk, vk, kplus, kminus, qplus, qminus = params
    a, b, c, d = amps
    vf, vkplus, vkminus, vqplus, vqminus = map(vel, [kf, kplus, kminus, qplus, qminus])
    A = vqminus*abs2(a)/vf
    B = vqplus*abs2(b)/vf

    if all
        if imag(kplus) > 0
            C = 0
        elseif imag(kplus) ≈ 0
            C = vkplus*abs2(c)*(abs2(uk) - abs2(vk))/vf
        else
            println("im kplus < 0")
        end
        if imag(kminus) < 0
            D = 0
        elseif imag(kminus) ≈ 0
            D = vkminus*abs2(d)*(abs2(uk) - abs2(vk))/vf
        else
            println("im kminus > 0")
        end
        return convert(Vector{Float64}, [A, B, C, D])
    else
        return convert(Vector{Float64}, [A, B])
    end
end

function calc_coeffs(store, Z, points, energies)
    H = ħ*vel(kf)*Z  # [eV*s*m/s = eV*m]
    for i in 1:points
        params = get_params(energies[i])
        amps = calc_amplitudes(energies[i], H, params)
        store[i,:] = coeffs(amps, params)
    end
    return store
end

function fermi(E, V, T)
    # V : [eV/e]
    return 1/(1 + exp((E - V)/(kb*T)))
end

function createintegrand(V, T, H)
    function integrand(E)
        params = get_params(E)
        amps = calc_amplitudes(E, H, params)
        A, B = coeffs(amps, params, false)
        return (fermi(E, V, T) - fermi(E, 0, T))*(1 + A - B)
    end
    return integrand
end

function calc_current(Z, T, V)
    # without prefactor
    H = ħ*vel(kf)*Z  # [eV*s*m/s = eV*m]
    fun = createintegrand(V, T, H)
    integral, err = quadgk(fun, -3*kb*T, V+3*kb*T, rtol=1e-5)
    return integral
end

function initiateglobals(Tc)
    global kb = 8.62e-5  # [eV/K]
    global Δ = 1.76*kb*Tc  # [eV]
    global mc2 = 0.511e6  # [eV]
    global ħ = 6.58e-16  # [eV*s]
    global c = 3e8  # [m/s]
    global μ = 11.7  # [eV]
    global kf = √(2*mc2*μ/(ħ*c)^2)
end

function plotcoeffs()
    Tc = 1.2  # [K]
    initiateglobals(Tc)
    points = 500
    energies = range(1e-4Δ, 5*Δ, points)
    Zvals = [0, 0.3, 1, 3] # this seems to work...?
    coeffs = zeros(Float64, points, 4, length(Zvals))
    for (i, Z) in enumerate(Zvals)
        coeffs[:,:,i] = calc_coeffs(coeffs[:,:,i], Z, points, energies)
    end
    plots = collect(plot(energies./Δ, collect(coeffs[:, i, j] for i in 1:4), dpi=300, title=L"Z="*string(Zvals[j])) for j in 1:length(Zvals))
    display(plot(plots..., layout=(2,2), label=[L"A" L"B" L"C" L"D"], xlabel=L"E/\Delta"))
    # savefig("coeffs.png")
end

function plotcurrent()
    Tc = 1.2  # [K]
    initiateglobals(Tc)
    Zvals = [0, 0.5, 1, 50] # this seems to work...?
    plotlabels = reshape([L"Z = "*string(Z) for Z in Zvals], (1,length(Zvals)))
    T = 50e-3
    points = 300
    voltages = range(1e-5*Δ, 3*Δ, points)
    currents = [[calc_current(Z, T, V).*(1+Z^2) for V in voltages] for Z in Zvals]
    h = voltages[2] - voltages[1]
    cond = [(currents[i][3:end] - currents[i][1:end-2])./(2*h) for i in 1:length(Zvals)]
    currentplot = plot(voltages./Δ, currents./Δ, ylabel=L"eIR_N/\Delta", dpi=300)
    condplot = plot(voltages[2:end-1]./Δ, cond, ylabel=L"R_N\frac{dI}{dV}", dpi=300)
    display(plot(currentplot, condplot, xlabel=L"eV/\Delta", labels=plotlabels))
    # savefig("current.png")
end

function lab()
    Tc = 1.2  # [K]
    initiateglobals(Tc)
    points = 100
    Z = 10
    mid = 430
    stop = 750
    len = stop - mid + 1
    temps = [0.05, 0.1, 0.25, 0.35, 0.5, 0.75, 1]
    index = [1, 3, 5, 7]
    n = length(index)
    currlab = zeros(len, n)
    currbtk = zeros(points, n)
    biases = zeros(len, n)
    voltages = range(1e-5*Δ, 3.5*Δ, points)
    for i in 1:n
        data = matread("00"*string(index[i])*".mat")
        currlab[:, i] = data["current"][1, mid:stop]
        currbtk[:, i] = [calc_current(Z, temps[index[i]], V) for V in voltages]
        biases[:, i] = data["bias"]["voltArray"][1, mid:stop]
    end
    biases .-= biases[1, 1]
    labels = reshape(["T = $t K" for t in temps[index]], (1, n))
    p1 = plot(biases, currlab, xaxis = "V", yaxis="I", lw=2, labels=labels, dpi=300)
    p1 = plot!(twinx(), voltages, currbtk, ls=:dash, lw=2, legend=false, showaxis=false, dpi=300)
    condlab = [(currlab[3:end, i] - currlab[1:end-2, i])./(2*(biases[2, i]-biases[1, i])) for i in 1:n]
    condbtk = [(currbtk[3:end, i] - currbtk[1:end-2, i])./(2*(voltages[2]-voltages[1])) for i in 1:n]
    p2 = plot(biases[2:end-1, :], condlab, xaxis = "V", yaxis=L"\frac{dI}{dV}", lw=2, labels=labels, dpi=300)
    p2 = plot!(twinx(), voltages[2:end-1], condbtk, ls=:dash, lw=2, legend=false, showaxis=false, dpi=300)
    display(plot(p1, p2, layout=(2,1)))
    #savefig("lab.png")
end
end
