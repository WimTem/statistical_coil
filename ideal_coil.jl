using Plots, Statistics, Random

function generate_coil(n=100)
    X, Y = 0, 0
    ang = 0
    A = zeros(n,2)
    for i = 1:n
        r = bitrand(1)
        if r[1] == 0
            ang += pi/3
        else
            ang -= pi/3
        end
        X += cos(ang)
        Y += sin(ang)
        A[i,:] = [X, Y]
    end
    return vcat([0 0], A)
end

n = 500
coil = generate_coil(n)
X = sum(coil[:,1])/n
Y = sum(coil[:,2])/n
plot(coil[:,1], coil[:,2], axis=nothing, border=:none, aspect_ratio=:equal, label="Coil", legend=:outertopleft)
plot!([0], [0], markershape=:hexagon, markersize=3, label="Start")
plot!([coil[end,1]], [coil[end,2]], markershape=:hexagon, markersize=3, label="End")
plot!([X], [Y], markershape=:hexagon, markersize=5, label="CoM")
savefig("coil.png")
function coil_radius(x)
    return distance([0,0], x[end,:])
end

function coil_center_of_mass(x)
    return sum(x,dims=1)/size(x)[1]
end

function distance(x, y)
    return ((y[2] - x[2])^2 + (y[1] - x[1])^2) |> sqrt
end

function coil_gyration_radius(x)
    n = size(x)[1]
    M = coil_center_of_mass(x)
    d = zeros(n)
    for i = 1:n
        d[i] = distance(x[i,:], M)
    end
    return sum(d)/n
end


# Theory: r^2 = 6.s^2
# 1000 runs with coil length = 500
A = zeros(1000, 3)
for i = 1:1000
    coil = generate_coil(500)
    r = coil_radius(coil)
    s = coil_gyration_radius(coil)
    A[i,1:2] = [r, s]
    A[i, 3] = 6*s^2
end

err = A[:,3] .- A[:,1].^2

mu1 = Statistics.mean(err)
std1 = Statistics.std(err)
result = []

#Remove outliers
for i = 1:length(err)
    if err[i] >= mu1-3*std1
        if err[i] <= mu1+3*std1
            x = (err[i]-mu1)/std1
            push!(result, x)
        end
    end
end
println(size(err))
println(size(result))

p1 = scatter(A[:,1].^2, label="r^2", yaxis=:log)
p2 = scatter!(A[:,3], label="6.s^2", yaxis=:log)
p3 = histogram(result, normalize=true)
savefig("results.png")
plot(p2, p3, layout=(2,1))
