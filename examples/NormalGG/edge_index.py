import numpy as np
from scipy.sparse import csr_matrix
import pandas as pd
import h5py

def generate_edge_index(rve_size, grain_distribution, include_self_loops=False):
    """
    生成 edge_index 表示的图结构。
    """
    nx, ny, nz = rve_size
    gx, gy, gz = grain_distribution
    grain_size = (nx // gx, ny // gy, nz // gz)
    
    # 构建晶粒分布
    grain_map = np.zeros((nx, ny, nz), dtype=int)
    grain_id = 0
    for i in range(0, nx, grain_size[0]):
        for j in range(0, ny, grain_size[1]):
            for k in range(0, nz, grain_size[2]):
                grain_map[i:i + grain_size[0], j:j + grain_size[1], k:k + grain_size[2]] = grain_id
                grain_id += 1

    # 构建邻接关系（edge_index）
    node_count = gx * gy * gz
    adjacency_list = []

    def get_linear_index(x, y, z):
        return x * gy * gz + y * gz + z

    for x in range(gx):
        for y in range(gy):
            for z in range(gz):
                current_node = get_linear_index(x, y, z)
                neighbors = []

                # 检查 6-邻域
                if x > 0:  # 左
                    neighbors.append(get_linear_index(x - 1, y, z))
                if x < gx - 1:  # 右
                    neighbors.append(get_linear_index(x + 1, y, z))
                if y > 0:  # 前
                    neighbors.append(get_linear_index(x, y - 1, z))
                if y < gy - 1:  # 后
                    neighbors.append(get_linear_index(x, y + 1, z))
                if z > 0:  # 下
                    neighbors.append(get_linear_index(x, y, z - 1))
                if z < gz - 1:  # 上
                    neighbors.append(get_linear_index(x, y, z + 1))

                for neighbor in neighbors:
                    adjacency_list.append((current_node, neighbor))

                # 添加自环
                if include_self_loops:
                    adjacency_list.append((current_node, current_node))

    # 转换为稀疏矩阵形式的边列表
    rows, cols = zip(*adjacency_list)
    edge_index = (np.array(rows, dtype=int), np.array(cols, dtype=int))
    return edge_index


def save_edge_index_to_hdf5(edge_index, filename):
    """
    将 edge_index 保存为 HDF5 文件。
    """
    with h5py.File(filename, "w") as f:
        f.create_dataset("edge_index/row", data=edge_index[0], compression="gzip")
        f.create_dataset("edge_index/col", data=edge_index[1], compression="gzip")

# 参数设置
rve_size = (30, 30, 30)  # RVE 的单元尺寸
grain_distribution = (15, 15, 15)  # 晶粒分布
include_self_loops = True  # 是否包含自环

# 生成 edge_index
edge_index = generate_edge_index(rve_size, grain_distribution, include_self_loops)

# 保存到 CSV 文件和 HDF5 文件
hdf5_filename = "edge_index.hdf5"
save_edge_index_to_hdf5(edge_index, hdf5_filename)

# 打印结果
print(f"Edge index 已保存为 HDF5 文件：{hdf5_filename}")
