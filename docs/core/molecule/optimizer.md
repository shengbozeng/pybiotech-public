# Optimizer — RDKit 3D 构象生成 / 力场优化 / 高保真二进制序列化（IPC 友好）

`Optimizer` 将**构象生成 + 力场优化**与**分子对象的高保真二进制传输**合在一起.
适用于协同计算、跨进程/跨应用传递,RDKit `Mol` 对象（如：Pipe/Socket/共享内存/文件）。

- 3D：ETKDGv3 嵌入 → 优先 **MMFF94**，失败回退 **UFF**
- 序列化：基于 **RDKit 原生二进制**（`MolToBinary`/`MolFromBinary`），**完整保留全部构象、坐标、手性与（可选）属性**
- 容器协议（v2）：批量打包为一个 blob，含 **版本/标志/长度前缀/CRC32**，可校验、可靠、向后兼容

> 典型场景：跨进程（Pipe/Socket/共享内存/文件）传递 `Mol`，不丢任何结构/构象/属性；另一端直接还原使用。

---

## 目录

- [Optimizer — RDKit 3D 构象生成 / 力场优化 / 高保真二进制序列化（IPC 友好）](#optimizer--rdkit-3d-构象生成--力场优化--高保真二进制序列化ipc-友好)
  - [目录](#目录)
  - [安装与依赖](#安装与依赖)
  - [快速上手](#快速上手)
  - [API 参考](#api-参考)
    - [`Optimizer.mmff_optimize`](#optimizermmff_optimize)
    - [`Optimizer.uff_optimize`](#optimizeruff_optimize)
    - [`Optimizer.embed_and_optimize`](#optimizerembed_and_optimize)
    - [`Optimizer.to_serialize`](#optimizerto_serialize)
    - [`Optimizer.to_unserialize(serialized_mols: bytes) -> List[Optional[Mol]]`](#optimizerto_unserializeserialized_mols-bytes---listoptionalmol)
  - [容器协议说明（v2，兼容 v1）](#容器协议说明v2兼容-v1)
  - [使用建议与最佳实践](#使用建议与最佳实践)
  - [错误处理与边界情况](#错误处理与边界情况)
  - [性能与兼容性提示](#性能与兼容性提示)
  - [示例：文件/Socket/共享内存传输](#示例文件socket共享内存传输)
  - [常见问题](#常见问题)

---

## 安装与依赖

- Python 3.8+
- RDKit（建议 2022.03+；已在 2024.09.x 验证）
- 标准库：`io`、`struct`、`zlib`、`typing`

> 共享内存示例需要 Python 3.8+ 的 `multiprocessing.shared_memory`。

---

## 快速上手

```python
from rdkit import Chem
from optimizer import Optimizer  # 假设文件名为 optimizer.py

# 1) 构造分子
mol = Chem.MolFromSmiles("CCO")

# 2) 生成 3D + 力场优化
mol3d, ok, meta = Optimizer.embed_and_optimize(
    mol,
    random_seed=123,
    max_ff_iters=800
)
print("收敛:", ok, "| 方法:", meta.get("method"), "| 能量:", meta.get("energy"))

# （可选）把优化信息写入属性（便于传输后审计）
if meta.get("method") is not None:
    mol3d.SetProp("opt.method", meta["method"])
if meta.get("energy") is not None:
    mol3d.SetProp("opt.energy", str(meta["energy"]))
mol3d.SetProp("opt.ok", "1" if ok else "0")

# 3) 批量序列化（Mol 或 None 均可；含属性、带 CRC 校验）
blob = Optimizer.to_serialize([mol3d, None], include_props=True, with_checksum=True)

# 4) 反序列化（顺序/占位一一对应）
mols = Optimizer.to_unserialize(blob)
print(type(mols[0]).__name__, mols[1])  # Mol None
```

---

## API 参考

### `Optimizer.mmff_optimize`

(m: Mol, max_ff_iters: int = 500) -> Optional[Dict[str, Any]]` （optional）

**作用**

- 使用 **MMFF94** 对分子进行力场最小化；不支持/失败返回 `None`。

**参数**

- `m`：`rdkit.Chem.Mol`（通常需已含 3D 坐标与显式氢）
- `max_ff_iters`：最大迭代次数（默认 500）

**返回（成功）**

```python
{
  "method": "MMFF94",
  "steps": <Minimize 返回码>,  # 0=收敛；>0=未收敛（注意：不是“实际步数”）
  "energy": <最终能量 float>
}
```

**示例**

```python
res = Optimizer.mmff_optimize(mol, max_ff_iters=800)
if res:
    print(res["method"], res["energy"], res["steps"])
```

---

### `Optimizer.uff_optimize`

(m: Mol, max_ff_iters: int = 500) -> Optional[Dict[str, Any]]` （optional）

**作用**

- 使用 **UFF** 进行最小化；异常或不可用时返回 `None`。

**参数/返回** 与 `mmff_optimize` 类似，其中 `method` 为 `"UFF"`。

**示例**

```python
res = Optimizer.uff_optimize(mol, max_ff_iters=1000)
if res:
    print(res["method"], res["energy"], res["steps"])
```

---

### `Optimizer.embed_and_optimize`

(mol,
*,
max_embed_attempts=1000,
random_seed=0xC0FFEE,
use_small_ring_torsions=True,
use_macrocycle_torsions=True,
prune_rms_thresh=0.1,
max_ff_iters=500)
->
(Mol, bool, Dict[str, Any])

**作用**

- 使用 **ETKDGv3** 生成初始 3D 构象，优先 **MMFF94** 优化，失败回退 **UFF**。

**流程**

- 复制入参（不修改原对象）
- `SanitizeMol`、`AssignStereochemistry`
- `AddHs(addCoords=True)`
- `ETKDGv3` 嵌入（清空旧构象，仅保留 1 个新构象）
- `MMFF94` → 不支持则 `UFF`

**关键参数**

- `random_seed`：固定可复现，`-1` 为完全随机
- `prune_rms_thresh`：去冗阈值（Å）
- `max_ff_iters`：力场最小化最大迭代

**返回**

- `optimized_mol: Mol`：含显式氢与 3D 构象；失败也返回当前工作副本便于诊断
- `ok: bool`：是否达到收敛（`Minimize` 返回 0）
- `meta: Dict[str, Any]`：诊断信息，例如

```python
{
  "stage": "optimize",   # or "init"/"sanitize"/"embed"
  "method": "MMFF94",    # or "UFF"
  "energy": 12.3456,
  "steps": 0,            # 0=收敛；>0=未收敛（返回码）
  "message": "..."
}
```

**示例**

```python
mol3d, ok, meta = Optimizer.embed_and_optimize(
    mol,
    random_seed=42,
    max_embed_attempts=2000,
    max_ff_iters=1000
)
```

---

### `Optimizer.to_serialize`

(
  mols: Iterable[Optional[Mol]],
  *,
  include_props: bool = True,
  with_checksum: bool = True
  )
->
bytes

**作用**

- 将**一组** `Mol`（或 `None` 占位）打包为**单一二进制容器**（v2），适合跨进程/跨应用传输或持久化。

**特性**

- **一一对应**：输入的每个元素输出为一条记录；`None` 会写为空记录（占位）
- **保真**：RDKit 原生二进制，完整保留所有构象/坐标/手性；`include_props=True` 时保留属性
- **强校验**：`with_checksum=True` 时为**每条记录**添加 CRC32 校验
- **版本记录**：容器头部写入 `rdBase.rdkitVersion`

**返回**

- `bytes`：容器 blob

**示例**

```python
blob = Optimizer.to_serialize(
    [mol_a, None, mol_b],
    include_props=True,
    with_checksum=True
)
```

---

### `Optimizer.to_unserialize(serialized_mols: bytes) -> List[Optional[Mol]]`

**作用**

- 还原容器 blob 为 `List[Optional[Mol]]`（**顺序/占位一一对应**）。

**兼容性**

- v2 容器：`b"RDKB\x02"`（含版本/flags/CRC）
- v1 容器：`b"RDKB\x01"`（无版本/flags/CRC）
- 若未匹配容器头：尝试将整个数据当作**裸 RDKit 单体二进制**；成功则返回长度为 1 的列表

**示例**

```python
mols = Optimizer.to_unserialize(blob)
for i, m in enumerate(mols):
    print(i, type(m).__name__ if m else None)
```

---

## 容器协议说明（v2，兼容 v1）

```text
字节序：大端（network order）

v2 Header
---------
5B  magic     = "RDKB\x02"
2B  vlen_be   = RDKit 版本字符串长度（uint16）
vlen  version = RDKit 版本 UTF-8（示例："2024.09.6"）
1B  flags     = bit0: include_props, bit1: with_checksum，其余保留为 0
4B  count_be  = 记录数 N（uint32）

v2 Records × N
--------------
4B  len_be    = 记录长度（int32）
                -1 → None 占位
                >=0 → 后续 payload 长度
len  payload  = RDKit Mol 二进制（当 len>0）
4B  crc_be    = （可选）CRC32（uint32；仅当 flags.with_checksum=1）
                None 记录 CRC 记作 0

v1 Header（兼容读取）
-------------------
5B  magic     = "RDKB\x01"
4B  count_be  = 记录数 N
records      = [len_be|payload] × N （无 CRC）

裸 RDKit 单体（兼容读取）
-----------------------
无法识别容器头时，尝试将整个 blob 作为单体 RDKit 二进制解析。
```

---

## 使用建议与最佳实践

- **复现性**：需要可复现实验时固定 `random_seed`；需要多样性时设为 `-1`
- **困难分子**：增大 `max_embed_attempts` 与 `max_ff_iters`；对大环保持 `use_macrocycle_torsions=True`
- **质量审计**：将 `meta` 写入属性（如 `opt.method/opt.energy/opt.ok`）并设置 `include_props=True`
- **批量/管道**：`to_serialize([...])` → 通过 Pipe/Socket/共享内存/文件发送 → `to_unserialize(...)`
- **校验**：跨进程传输建议保持 `with_checksum=True`（默认）
- **长期归档**：内部传输二进制最快最全；若需多年跨版本可读性，额外保留 SDF/JSON 备份

---

## 错误处理与边界情况

- `embed_and_optimize`

  - `mol is None` → `ValueError`
  - 嵌入失败或力场不可用 → 返回 `(work_copy, False, meta)`（不抛异常，`meta["stage"]/["message"]` 说明原因）
  - `meta["steps"]` 是**返回码**（0=收敛，>0=未收敛），不是迭代步数
- `to_serialize`

  - 输入序列为空 → `ValueError`
- `to_unserialize`

  - 空输入、截断、长度异常、CRC 不匹配 → `ValueError`
  - 兼容 v1 容器与“裸二进制单体”兜底

---

## 性能与兼容性提示

- **性能**：二进制序列化/反序列化通常显著快于文本格式，体积也更小；CRC32 为 C 实现，开销低
- **兼容性**：RDKit 二进制用于在线/IPC 非常稳妥；若面向超长期归档或跨很多 RDKit 版本，建议同时保留文本/JSON 备份
- **属性**：`include_props=True` 时，Mol/Atom/Bond 属性会被保留（字符串/数值建议转字符串或用 `SetDoubleProp`）

---

## 示例：文件/Socket/共享内存传输

**写入文件 → 读取**

```python
# 写
blob = Optimizer.to_serialize([mol_a, None, mol_b], include_props=True, with_checksum=True)
with open("mols.bin", "wb") as f:
    f.write(blob)

# 读
with open("mols.bin", "rb") as f:
    data = f.read()
mols = Optimizer.to_unserialize(data)
```

**通过 Socket 发送（核心思路：先发长度，再发数据）**

```python
# sender
blob = Optimizer.to_serialize(mols, include_props=True, with_checksum=True)
sock.sendall(len(blob).to_bytes(8, "big"))
sock.sendall(blob)

# receiver
import struct
def recv_exact(s, n):
    buf = bytearray()
    while len(buf) < n:
        chunk = s.recv(n - len(buf))
        if not chunk:
            raise ConnectionError("socket closed")
        buf.extend(chunk)
    return bytes(buf)

hdr = recv_exact(sock, 8)
length = struct.unpack(">Q", hdr)[0]
payload = recv_exact(sock, length)
mols = Optimizer.to_unserialize(payload)
```

**共享内存（摘要）**

```python
from multiprocessing import shared_memory
blob = Optimizer.to_serialize(mols, include_props=True, with_checksum=True)

# producer
shm = shared_memory.SharedMemory(create=True, size=len(blob))
shm.buf[:len(blob)] = blob
# 把 shm.name 和 len(blob) 通过管道/消息发送给消费者
# 消费完成后由 producer 做 shm.unlink()

# consumer
shm = shared_memory.SharedMemory(name=shm_name, create=False)
view = bytes(shm.buf[:blob_len])  # 或 memoryview 再转 bytes
mols = Optimizer.to_unserialize(view)
```

---

## 常见问题

**Q：二进制会丢“力场方法/能量/是否收敛”吗？**
A：这些并非 Mol 的内建字段（无论文本/二进制都不会自动保存）。如需保留，请在优化后写入 Mol 属性，并在序列化时使用 `include_props=True`。

**Q：为什么 `meta["steps"]` 不是实际步数？**
A：RDKit 的 `ForceField.Minimize()` 返回的是**返回码**：0=收敛，其它值=未收敛。实际迭代步数不直接暴露。

**Q：容器里能放很多分子吗？**
A：可以。建议按业务需要切批（例如 1k–10k/批），方便流式处理与失败重试。

**Q：CRC 校验能关吗？**
A：可通过 `with_checksum=False` 关闭以追求极致吞吐；跨进程/跨主机建议保持开启（默认）。
