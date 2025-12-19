import h5py
import argparse
from pathlib import Path

def list_available_timesteps(h5_file: str):
    with h5py.File(h5_file, 'r') as f:
        grp = '/CheckPoints/GrainConnections'
        if grp not in f:
            raise ValueError(f"Group {grp} not found in HDF5 file: {h5_file}")
        return sorted(int(k) for k in f[grp].keys())

def read_grain_connections(h5_file, timestep):
    # 打开HDF5文件
    with h5py.File(h5_file, 'r') as f:
        # 构造路径
        dataset_path = f"/CheckPoints/GrainConnections/{timestep}"
        
        # 检查数据集是否存在
        if dataset_path not in f:
            raise ValueError(f"Dataset {dataset_path} not found in HDF5 file.")
        
        # 读取数据
        data = f[dataset_path][:]
    
    # 解析扁平化数据
    grain_connections = {}
    i = 0
    while i < len(data):
        grain_id = int(data[i])  # 当前晶粒ID
        neighbor_count = int(data[i + 1])  # 邻居数量
        neighbors = [int(data[i + 2 + j]) for j in range(neighbor_count)]  # 邻居ID列表
        grain_connections[grain_id] = neighbors
        i += 2 + neighbor_count  # 跳到下一个晶粒的数据
    
    return grain_connections

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Read GrainConnections from HDF5 and print adjacency list")
    parser.add_argument("--file", "-f", default="NormalGG_output.h5", help="Path to HDF5 file (default: NormalGG_output.h5)")
    parser.add_argument("--timestep", "-t", type=int, default=None, help="Timestep to read (default: latest available)")
    args = parser.parse_args()

    # Resolve file path relative to this script's directory if a bare name is given
    h5_file = args.file
    if not Path(h5_file).exists():
        candidate = Path(__file__).resolve().parent / h5_file
        if candidate.exists():
            h5_file = str(candidate)

    # Determine timestep
    available = list_available_timesteps(h5_file)
    if not available:
        raise RuntimeError("No timesteps found under /CheckPoints/GrainConnections")
    if args.timestep is None:
        timestep = available[-1]
        print(f"[Info] Using latest timestep: {timestep}")
    else:
        timestep = args.timestep
        if timestep not in available:
            raise ValueError(f"Timestep {timestep} not found. Available: {available[:5]} ... {available[-5:]}")

    connections = read_grain_connections(h5_file, timestep)
    print(f"[Info] Read {len(connections)} grains at timestep {timestep} from {h5_file}")
    # 打印每个晶粒的连接关系
    for grain, neighbors in connections.items():
        print(f"Grain {grain} is connected to: {neighbors}")