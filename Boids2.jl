# Boids2.jl: メトリック型とトポロジカル型の群れモデル比較(修正版)
include("params2.jl")
include("obstacles2.jl")

using Random, LinearAlgebra, Statistics, Plots, Printf
gr()

# 再現性確保
Random.seed!(365)

# --- データ構造 ---
mutable struct Boid
    pos::Vector{Float64}
    vel::Vector{Float64}
end

# --- 初期化 ---
function initialize_boids()
    [Boid(rand(2) .* [WIDTH, HEIGHT], (rand(2) .- 0.5) .* max_speed) for _ in 1:N]
end

# --- トーラス境界（周期境界）計算 ---
function torus_vector(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    dx = b[1] - a[1]
    if dx > WIDTH/2; dx -= WIDTH; elseif dx < -WIDTH/2; dx += WIDTH; end

    dy = b[2] - a[2]
    if dy > HEIGHT/2; dy -= HEIGHT; elseif dy < -HEIGHT/2; dy += HEIGHT; end

    return [dx, dy]
end

function torus_distance(a::AbstractVector{<:Real}, b::AbstractVector{<:Real})
    return norm(torus_vector(a, b))
end

# --- 近傍個体の抽出 ---
# メトリック型
function neighbors_metric(i::Int, boids::Vector{Boid}, r::Float64)
    res = Int[]
    pi = boids[i].pos
    for j in 1:length(boids)
        if i == j; continue; end
        if torus_distance(pi, boids[j].pos) < r; push!(res, j); end
    end
    return res
end

# トポロジカル型
function neighbors_topological(i::Int, boids::Vector{Boid}, k::Int)
    pi = boids[i].pos
    dists = [(j, torus_distance(pi, boids[j].pos)) for j in 1:length(boids) if j != i]
    sorted = sort(dists, by = x -> x[2])
    m = min(k, length(sorted))
    return [sorted[t][1] for t in 1:m]
end

# --- 個体状態の更新 ---
function update_boids!(boids::Vector{Boid}, neighbor_fn::Function, neighbor_param, obstacles)
    Nloc = length(boids)
    new_vels = Vector{Vector{Float64}}(undef, Nloc)

    for i in 1:Nloc
        b = boids[i]
        neighbors_idx = neighbor_fn(i, boids, neighbor_param)
        cohesion, alignment = zeros(2), zeros(2)

        if !isempty(neighbors_idx)
            # 結合(Cohesion)と整列(Alignment)
            deltas = [torus_vector(b.pos, boids[j].pos) for j in neighbors_idx]
            mean_delta = sum(deltas) ./ length(deltas)
            avg_pos = b.pos .+ mean_delta
            avg_vel = sum(boids[j].vel for j in neighbors_idx) ./ length(neighbors_idx)

            cohesion .= (avg_pos - b.pos) .* w_cohesion
            alignment .= (avg_vel - b.vel) .* w_alignment
        end

        # 分離(Separation)
        close_idx = [j for j in 1:Nloc if i != j && torus_distance(b.pos, boids[j].pos) < separation_radius]
        separation = zeros(2)
        if !isempty(close_idx)
            sep_sum = zeros(2)
            for j in close_idx
                diff = -torus_vector(b.pos, boids[j].pos)
                sep_sum .+= diff ./ (norm(diff)^2 + 1e-5)
            end
            separation .= sep_sum .* w_separation
        end

        # 最終速度の算出（障害物回避含む）
        obstacle_force = avoid_obstacles(b, obstacles, avoidance_radius) * w_obstacle #障害物回避の重み
        new_vel = b.vel .+ cohesion .+ separation .+ alignment .+ obstacle_force

        vnorm = norm(new_vel)
        if vnorm > 1e-8; new_vel .+= (new_vel / vnorm) * propulsion_strength; end

        vnorm = norm(new_vel)
        if vnorm > max_speed; new_vel = new_vel ./ vnorm .* max_speed; end

        new_vels[i] = new_vel
    end

    # 位置の同期更新
    for i in 1:length(boids)
        b = boids[i]
        b.vel .= new_vels[i]
        b.pos .+= b.vel .* dt
        for d in 1:2
            limit = d == 1 ? WIDTH : HEIGHT
            if b.pos[d] < 0; b.pos[d] += limit; elseif b.pos[d] > limit; b.pos[d] -= limit; end
        end
    end
end

# --- 障害物回避 ---
function avoid_obstacles(b::Boid, obstacles, avoidance_radius)
    force = zeros(2)
    for ob in obstacles
        d = torus_vector(b.pos, ob.pos)
        dist = norm(d)
        if dist < ob.r + avoidance_radius; force .+= -d / (dist^2 + 1e-5); end
    end
    return force
end

# --- 可視化関数 ---
function draw_obstacles!(p, obstacles)
    for ob in obstacles
        col = ob.movable ? :red : :green
        scatter!(p, [ob.pos[1]], [ob.pos[2]], markersize = ob.r * 10,
                 markercolor = col, markerstrokecolor = :black, markerstrokewidth = 1.5, alpha = 0.8)
    end
end

function plot_pair(boids_metric, boids_topo, obstacles, t)
    # Metric Plot
    p1 = quiver([b.pos[1] for b in boids_metric], [b.pos[2] for b in boids_metric],
                quiver=([b.vel[1] for b in boids_metric], [b.vel[2] for b in boids_metric]),
                xlims=(0, WIDTH), ylims=(0, HEIGHT), aspect_ratio=1,
                title=@sprintf("Metric (r=%.1f) Step: %d", metric_radius, t), legend=false, color=:blue)
    draw_obstacles!(p1, obstacles)

    # Topological Plot
    p2 = quiver([b.pos[1] for b in boids_topo], [b.pos[2] for b in boids_topo],
                quiver=([b.vel[1] for b in boids_topo], [b.vel[2] for b in boids_topo]),
                xlims=(0, WIDTH), ylims=(0, HEIGHT), aspect_ratio=1,
                title=@sprintf("Topological (k=%d) Step: %d", k_topo, t), legend=false, color=:red)
    draw_obstacles!(p2, obstacles)

    return plot(p1, p2, layout=(1, 2), size=(1000, 450))
end

# --- 指標計算 ---
function calculate_order(boids::Vector{Boid})
    normalized_vels = [b.vel / norm(b.vel) for b in boids if norm(b.vel) > 1e-8]
    if isempty(normalized_vels); return 0.0; end
    return norm(sum(normalized_vels) / length(normalized_vels))
end

function mean_nn_distance(positions::AbstractMatrix)
    Nloc, d = size(positions)
    nn = zeros(Nloc)
    for i in 1:Nloc
        min_d = Inf
        pi = view(positions, i, :)
        for j in 1:Nloc
            if i == j; continue; end
            dj = norm(torus_vector(collect(pi), collect(view(positions, j, :))))
            if dj < min_d; min_d = dj; end
        end
        nn[i] = min_d
    end
    return mean(nn)
end

function cluster_size(boids::Vector{Boid}, R_cluster::Float64)
    Nloc = length(boids)
    visited = falses(Nloc)
    max_cluster = 0

    for i in 1:Nloc
        if visited[i]; continue; end

        stack = [i]
        visited[i] = true
        size = 1

        while !isempty(stack)
            v = pop!(stack)
            for j in 1:Nloc
                if visited[j]; continue; end
                if torus_distance(boids[v].pos, boids[j].pos) < R_cluster
                    push!(stack, j)
                    visited[j] = true
                    size += 1
                end
            end
        end

        max_cluster = max(max_cluster, size)
    end

    return max_cluster
end

function speed_stats(boids)
    speeds = [norm(b.vel) for b in boids]
    return mean(speeds), std(speeds)
end

# --- シミュレーション実行 ---
boids_init = initialize_boids()
boids_metric, boids_topo = deepcopy(boids_init), deepcopy(boids_init)

order_history_metric, order_history_topo = Float64[], Float64[]
nn_history_metric, nn_history_topo = Float64[], Float64[]
obstacles = Obstacle[]

cluster_history_metric = Int[]
cluster_history_topo = Int[]

mean_v_history_metric = Float64[]
std_v_history_metric  = Float64[]

mean_v_history_topo = Float64[]
std_v_history_topo  = Float64[]

# --- 障害物の追加 ---
#add_obstacle!(obstacles, create_static_obstacle([25.0, 25.0], 4.0))
add_obstacle!(obstacles, create_moving_obstacle([45.0, 25.0], 4.0, [-1.5, 0.0]))

# シミュレーションループ
anim = @animate for t in 1:T
    update_obstacles!(obstacles, dt, WIDTH, HEIGHT)
    update_boids!(boids_metric, neighbors_metric, metric_radius, obstacles)
    update_boids!(boids_topo, neighbors_topological, k_topo, obstacles)

    push!(order_history_metric, calculate_order(boids_metric))
    push!(order_history_topo, calculate_order(boids_topo))

    pos_metric = reduce(hcat, [b.pos for b in boids_metric])'
    pos_topo = reduce(hcat, [b.pos for b in boids_topo])'

    push!(nn_history_metric, mean_nn_distance(pos_metric))
    push!(nn_history_topo, mean_nn_distance(pos_topo))

    push!(cluster_history_metric,cluster_size(boids_metric, cluster_radius))
    push!(cluster_history_topo,cluster_size(boids_topo, cluster_radius))

    mean_v_m, std_v_m = speed_stats(boids_metric)
    mean_v_t, std_v_t = speed_stats(boids_topo)
    
    push!(mean_v_history_metric, mean_v_m)
    push!(std_v_history_metric,  std_v_m)

    push!(mean_v_history_topo, mean_v_t)
    push!(std_v_history_topo,  std_v_t)

    plot_pair(boids_metric, boids_topo, obstacles, t)
end

# 結果保存
gif(anim, "Boids-C.gif", fps=20)

p_order = plot(order_history_metric, label="Metric", xlabel="Step", ylabel="Order Parameter",
               title="Order Parameter Comparison", linewidth=2, color=:blue)
plot!(p_order, order_history_topo, label="Topological", color=:red, linewidth=2)
savefig(p_order, "order-C.svg")

p_nn = plot(nn_history_metric, label="Metric", xlabel="Step", ylabel="Mean NN Distance",
            title="Mean NN Distance Comparison", linewidth=2, color=:blue)
plot!(p_nn, nn_history_topo, label="Topological", color=:red, linewidth=2)
savefig(p_nn, "nn-C.svg")

p_cluster = plot(cluster_history_metric, label="Metric",
    xlabel="Step", ylabel="Max Cluster Size",
    title="Maximum Cluster Size Comparison", linewidth=2)

plot!(p_cluster, cluster_history_topo, label="Topological", linewidth=2)
savefig(p_cluster, "cluster-C.svg")

p_std = plot(std_v_history_metric, label="Metric", xlabel="Step", ylabel="Speed Std",
            title="Speed Dispersion", linewidth=2)
plot!(p_std, std_v_history_topo, label="Topological", linewidth=2)
savefig(p_std, "speed_std-C.svg")

println("Finished: result files saved.")