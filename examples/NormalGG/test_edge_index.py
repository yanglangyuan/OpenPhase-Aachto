#!/usr/bin/env python
"""
读取并验证HDF5文件中的EdgeIndex格式
"""
import h5py
import numpy as np
import argparse
from pathlib import Path

def list_available_timesteps(h5_file: str, dataset_name: str = "EdgeIndex"):
    """列出可用的时间步"""
    with h5py.File(h5_file, 'r') as f:
        grp = f'/CheckPoints/{dataset_name}'
        if grp not in f:
            raise ValueError(f"Group {grp} not found in HDF5 file: {h5_file}")
        return sorted(int(k) for k in f[grp].keys())

def read_edge_index(h5_file: str, timestep: int):
    """
    读取EdgeIndex (row, col) 格式
    新的结构: /CheckPoints/EdgeIndex/{timestep}/row 和 /CheckPoints/EdgeIndex/{timestep}/col
    返回：(row_array, col_array)
    """
    with h5py.File(h5_file, 'r') as f:
        row_path = f"/CheckPoints/EdgeIndex/{timestep}/row"
        col_path = f"/CheckPoints/EdgeIndex/{timestep}/col"
        
        if row_path not in f or col_path not in f:
            raise ValueError(f"EdgeIndex datasets not found at timestep {timestep}")
        
        row = f[row_path][:]
        col = f[col_path][:]
    
    return row.astype(int), col.astype(int)

def read_grain_connections(h5_file: str, timestep: int):
    """
    读取GrainConnections扁平化格式
    """
    with h5py.File(h5_file, 'r') as f:
        dataset_path = f"/CheckPoints/GrainConnections/{timestep}"
        if dataset_path not in f:
            raise ValueError(f"Dataset {dataset_path} not found in HDF5 file.")
        data = f[dataset_path][:]
    
    grain_connections = {}
    i = 0
    while i < len(data):
        grain_id = int(data[i])
        neighbor_count = int(data[i + 1])
        neighbors = [int(data[i + 2 + j]) for j in range(neighbor_count)]
        grain_connections[grain_id] = neighbors
        i += 2 + neighbor_count
    
    return grain_connections

def verify_edge_index_consistency(grain_connections: dict, edge_row: np.ndarray, edge_col: np.ndarray) -> bool:
    """
    验证EdgeIndex和GrainConnections的一致性
    对于无向图，应该有 (i,j) 当且仅当有 (j,i)
    """
    # 从edge_index构建邻接表
    edge_adjacency = {}
    for src, dst in zip(edge_row, edge_col):
        if src not in edge_adjacency:
            edge_adjacency[src] = []
        edge_adjacency[src].append(dst)
    
    # 检查一致性
    all_consistent = True
    for grain_id, neighbors in grain_connections.items():
        edge_neighbors = sorted(edge_adjacency.get(grain_id, []))
        grain_neighbors = sorted(neighbors)
        
        if edge_neighbors != grain_neighbors:
            print(f"❌ Grain {grain_id}:")
            print(f"   GrainConnections: {grain_neighbors}")
            print(f"   EdgeIndex: {edge_neighbors}")
            all_consistent = False
    
    return all_consistent

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read and verify EdgeIndex from HDF5")
    parser.add_argument("--file", "-f", default="NormalGG_output.h5", help="Path to HDF5 file")
    parser.add_argument("--timestep", "-t", type=int, default=None, help="Timestep to read (default: latest)")
    parser.add_argument("--verify", "-v", action="store_true", help="Verify consistency between EdgeIndex and GrainConnections")
    args = parser.parse_args()
    
    h5_file = args.file
    if not Path(h5_file).exists():
        candidate = Path(__file__).resolve().parent / h5_file
        if candidate.exists():
            h5_file = str(candidate)
    
    # 确定时间步
    available = list_available_timesteps(h5_file)
    if not available:
        raise RuntimeError("No timesteps found")
    
    if args.timestep is None:
        timestep = available[-1]
        print(f"[Info] Using latest timestep: {timestep}\n")
    else:
        timestep = args.timestep
        if timestep not in available:
            raise ValueError(f"Timestep {timestep} not found")
    
    # 读取EdgeIndex
    edge_row, edge_col = read_edge_index(h5_file, timestep)
    print(f"[EdgeIndex] Timestep {timestep}:")
    print(f"  - Row shape: {edge_row.shape}")
    print(f"  - Col shape: {edge_col.shape}")
    print(f"  - Number of edges: {len(edge_row)}")
    print(f"  - Node range: [{edge_row.min()}, {edge_row.max()}] (row), [{edge_col.min()}, {edge_col.max()}] (col)")
    print(f"\n  First 10 edges:")
    for i in range(min(10, len(edge_row))):
        print(f"    Edge {i}: {edge_row[i]} -> {edge_col[i]}")
    
    # 如果有验证选项，读取GrainConnections并验证一致性
    if args.verify:
        print(f"\n[Verification] Checking consistency with GrainConnections...")
        grain_connections = read_grain_connections(h5_file, timestep)
        if verify_edge_index_consistency(grain_connections, edge_row, edge_col):
            print("✅ All grains are consistent!")
        else:
            print("❌ Some inconsistencies found!")
