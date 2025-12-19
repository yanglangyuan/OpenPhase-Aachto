#!/usr/bin/env python
"""
示例：如何使用HDF5文件中的EdgeIndex数据进行机器学习/图神经网络
"""
import h5py
import numpy as np

def load_graph_data(h5_file: str, timestep: int):
    """
    从HDF5文件加载图数据，返回适用于PyTorch Geometric的格式
    """
    with h5py.File(h5_file, 'r') as f:
        # 读取边列表
        edge_row = f[f"/CheckPoints/EdgeIndex/{timestep}/row"][:]
        edge_col = f[f"/CheckPoints/EdgeIndex/{timestep}/col"][:]
        
        # 读取节点特征（晶粒体积和邻居数）
        grain_volumes = f[f"/CheckPoints/GrainVolumes/{timestep}"][:]
        grain_neighbors = f[f"/CheckPoints/GrainNeighbors/{timestep}"][:]
    
    # 构建edge_index (2, num_edges) 格式，适用于PyTorch Geometric
    edge_index = np.stack([edge_row, edge_col], axis=0)
    
    # 构建节点特征矩阵 (num_nodes, num_features)
    node_features = np.stack([grain_volumes, grain_neighbors], axis=1)
    
    return {
        'edge_index': edge_index,
        'node_features': node_features,
        'num_nodes': len(grain_volumes),
        'num_edges': len(edge_row)
    }

def example_pytorch_geometric_usage(h5_file: str, timestep: int):
    """
    示例：使用PyTorch Geometric加载数据
    """
    try:
        import torch
        from torch_geometric.data import Data
        
        data = load_graph_data(h5_file, timestep)
        
        # 转换为PyTorch Geometric的Data对象
        graph_data = Data(
            x=torch.tensor(data['node_features'], dtype=torch.float),
            edge_index=torch.tensor(data['edge_index'], dtype=torch.long)
        )
        
        print(f"\n[PyTorch Geometric] Graph at timestep {timestep}:")
        print(f"  - Number of nodes: {graph_data.num_nodes}")
        print(f"  - Number of edges: {graph_data.num_edges}")
        print(f"  - Node feature dimension: {graph_data.x.shape[1]}")
        print(f"  - Is undirected: {graph_data.is_undirected()}")
        print(f"\n  Node features (first 5 grains):")
        print(f"  [Volume, Neighbors]")
        for i in range(min(5, graph_data.num_nodes)):
            print(f"  Grain {i}: {graph_data.x[i].numpy()}")
        
        return graph_data
    except ImportError:
        print("PyTorch Geometric not installed. Install with: pip install torch-geometric")
        return None

def example_networkx_usage(h5_file: str, timestep: int):
    """
    示例：使用NetworkX加载数据
    """
    try:
        import networkx as nx
        
        data = load_graph_data(h5_file, timestep)
        
        # 创建NetworkX图
        G = nx.Graph()
        
        # 添加节点和特征
        for i in range(data['num_nodes']):
            G.add_node(i, 
                      volume=data['node_features'][i, 0],
                      neighbors=data['node_features'][i, 1])
        
        # 添加边
        edges = list(zip(data['edge_index'][0], data['edge_index'][1]))
        G.add_edges_from(edges)
        
        print(f"\n[NetworkX] Graph at timestep {timestep}:")
        print(f"  - Number of nodes: {G.number_of_nodes()}")
        print(f"  - Number of edges: {G.number_of_edges()}")
        print(f"  - Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")
        print(f"  - Is connected: {nx.is_connected(G)}")
        
        # 计算一些图统计量
        if nx.is_connected(G):
            print(f"  - Average clustering coefficient: {nx.average_clustering(G):.4f}")
            print(f"  - Graph density: {nx.density(G):.4f}")
        
        return G
    except ImportError:
        print("NetworkX not installed. Install with: pip install networkx")
        return None

def example_save_for_external_tools(h5_file: str, timestep: int, output_prefix: str = "grain_graph"):
    """
    示例：导出为其他工具可用的格式（CSV、edge list等）
    """
    data = load_graph_data(h5_file, timestep)
    
    # 保存边列表为CSV
    edge_list_file = f"{output_prefix}_t{timestep}_edges.csv"
    np.savetxt(edge_list_file, data['edge_index'].T, 
               delimiter=',', header='source,target', comments='', fmt='%d')
    print(f"\n[Export] Edge list saved to: {edge_list_file}")
    
    # 保存节点特征为CSV
    node_features_file = f"{output_prefix}_t{timestep}_nodes.csv"
    np.savetxt(node_features_file, data['node_features'], 
               delimiter=',', header='volume,neighbors', comments='', fmt='%.6e,%d')
    print(f"[Export] Node features saved to: {node_features_file}")
    
    return edge_list_file, node_features_file

if __name__ == "__main__":
    import argparse
    from pathlib import Path
    
    parser = argparse.ArgumentParser(description="Example ML/GNN usage of grain graph data")
    parser.add_argument("--file", "-f", default="NormalGG_output.h5", help="HDF5 file path")
    parser.add_argument("--timestep", "-t", type=int, default=0, help="Timestep to load")
    parser.add_argument("--export", action="store_true", help="Export to CSV files")
    args = parser.parse_args()
    
    h5_file = args.file
    if not Path(h5_file).exists():
        candidate = Path(__file__).resolve().parent / h5_file
        if candidate.exists():
            h5_file = str(candidate)
    
    print(f"Loading graph data from: {h5_file}")
    print(f"Timestep: {args.timestep}\n")
    
    # 基础数据加载
    data = load_graph_data(h5_file, args.timestep)
    print(f"[Basic Info]")
    print(f"  - Number of nodes (grains): {data['num_nodes']}")
    print(f"  - Number of edges (connections): {data['num_edges']}")
    print(f"  - Average degree: {data['num_edges'] / data['num_nodes']:.2f}")
    
    # NetworkX示例
    example_networkx_usage(h5_file, args.timestep)
    
    # PyTorch Geometric示例
    example_pytorch_geometric_usage(h5_file, args.timestep)
    
    # 导出选项
    if args.export:
        example_save_for_external_tools(h5_file, args.timestep)
