# =========================================================
# obstacles.jl: 障害物オブジェクトの定義と制御
# =========================================================

"""
Obstacle 構造体
- pos: 座標 [x, y]
- r:   半径（衝突判定に利用）
- vel: 速度ベクトル [vx, vy]
- movable: 移動の有無 (true/false)
"""
mutable struct Obstacle
    pos::Vector{Float64}
    r::Float64
    vel::Vector{Float64}
    movable::Bool
end

# --- 静的障害物の生成 ---
function create_static_obstacle(pos::Vector{Float64}, r::Float64)
    return Obstacle(pos, r, [0.0, 0.0], false)
end

# --- 動的障害物の生成 ---
function create_moving_obstacle(pos::Vector{Float64}, r::Float64, vel::Vector{Float64})
    return Obstacle(pos, r, vel, true)
end

# --- シミュレーションへの追加 ---
function add_obstacle!(obs_list::Vector{Obstacle}, obs::Obstacle)
    push!(obs_list, obs)
end

# --- 状態更新（移動と周期境界処理） ---
function update_obstacles!(obs_list::Vector{Obstacle}, dt::Float64, WIDTH::Float64, HEIGHT::Float64)
    for ob in obs_list
        if ob.movable
            ob.pos .+= ob.vel .* dt

            # 周期境界処理: mod関数を使用してマイナス速度時も安全にループさせる
            ob.pos[1] = mod(ob.pos[1], WIDTH)
            ob.pos[2] = mod(ob.pos[2], HEIGHT)
        end
    end
end