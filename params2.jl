# params2.jl

# --- シミュレーション全体パラメータ ---
# 個体数
const N = 100
# 表示幅
const WIDTH, HEIGHT = 50.0, 50.0
# 1ステップの時間
const dt = 0.1
# ステップ数
const T = 500
# 最大速度
const max_speed = 2.0
# 微小推進項
const propulsion_strength = 0.01

# --- 行動ルールの重み ---
# 結合
const w_cohesion = 0.01
# 整列
const w_alignment = 0.05
# 分離
const w_separation = 0.1

# --- 近傍条件 ---
# メトリックが作用する距離（半径）
const metric_radius = 10.0
# 分離が作用する距離（半径）
const separation_radius = 5.0
# トポロジカルで作用する個体数
const k_topo = 6

# --- 最大クラスタサイズでの条件 ---
# クラスタ判定距離
const cluster_radius = 4.0

# --- 障害物 ---
# 障害物回避作用半径
const avoidance_radius = 5.0
# 障害物回避の重み
const w_obstacle = 1.5