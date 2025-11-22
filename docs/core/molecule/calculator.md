# Molecule Calculator

`MolCaculator` 集中了一系列基于 RDKit/NumPy 的静态工具，用来获取化合物的几何、表面、指纹、环及药效团相关特征，适合围绕单个 `Chem.Mol` 在数据清洗或特征工程中的重复调用。✨

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `MolCaculator` | class | 仅包含静态方法的“分子计算器”，封装常见的几何、表面、指纹、药效团与电荷相关计算。 |
| `calc_surface_area` | static method | 通过 `rdFreeSASA` 计算总表面积，要求提前嵌入 3D。 |
| `calc_hydrophobic_surface_area` | static method | 计算疏水原子部分的 SASA，可选高精度原子感知。 |
| `calc_diameter` | static method | 计算 3D 坐标下分子的最大原子间距离。 |
| `calc_sssr` | static method | 获取 SSSR 环集、环尺寸分组与 atom/bond 环计数等结构化摘要。 |
| `calc_logp` | static method | 通过 `Crippen.MolLogP` 计算 LogP。 |
| `calc_morgan_fringerprint` | static method | 生成 Morgan 指纹并以 NumPy 数组形式返回。 |
| `calc_maccs_fingerprint` | static method | 生成 MACCS 指纹并转成 NumPy 数组。 |
| `calc_crippen_contribs` | static method | 计算每个原子的 Crippen 贡献（LogP 与 MolRefr）。 |
| `get_ring_atoms` | static method | 列出所有环原子的索引，可选 1-based 输出；支持 `ringinfo` 或 `GetSymmSSSR`。 |
| `get_hydrophobic_clusters` | static method | 依据 3D 距离和邻居元素谓词，聚类疏水碳原子。 |
| `approximate_oe_hydrophobe` | static method | 模拟 OEShape 的疏水簇识别，支持电荷/芳香性/计数等过滤。 |
| `get_anions_cations` | static method | 按 formal charge 拆分阴/阳离子原子索引。 |
| `get_hbond_acceptors` | static method | 通过 SMARTS 匹配氢键受体原子。 |
| `get_hbond_donors` | static method | 通过 SMARTS 匹配氢键供体原子。 |
| `calc_pharmacophore_features` | static method | 组合多个子模块构建药效团特征、计数与原子列表。 |
| `calc_intra_pharmacophore_distance` | static method | 计算同一药效团内部的距离矩阵。 |
| `calc_inter_pharmacophore_distance` | static method | 计算两个药效团之间的距离矩阵。 |

---

### `MolCaculator.calc_surface_area(mol: Chem.Mol) -> float`

**功能说明**

- 返回完整分子的表面积（SASA），使用 `rdFreeSASA.CalcSASA`，要求包含显式氢与有效 3D 坐标。

**参数**

- `mol`：含 3D 构象的 RDKit `Chem.Mol`；生产级调用前请执行 `Chem.AddHs` + Embed + Optimize。

**返回类型**

- `float`：单位为 Å² 的表面积总和。

**异常**

- 转换/计算异常会捕获并重新抛出 `Exception`。

**注意事项**

- ⚠️ 需预先调用 `mol.UpdatePropertyCache(strict=False)`、Embed、MMFF/UFF 优化等；否则 `rdFreeSASA` 会因缺少坐标而失败。

**示例**

```python
from rdkit import Chem
from pybiotech.core.molecule.calculator import MolCaculator

mol = Chem.MolFromSmiles("CCO")
mol = Chem.AddHs(mol)
Chem.AllChem.EmbedMolecule(mol)
area = MolCaculator.calc_surface_area(mol)
```

---

### `MolCaculator.calc_hydrophobic_surface_area(mol: Chem.Mol, is_high_precise: bool = False) -> float`

**功能说明**

- 计算指定为疏水的原子所占的 SASA；可选 `is_high_precise=True` 以扩展匹配至卤素/杂原子。

**参数**

- `mol`：已嵌入 3D 的分子。
- `is_high_precise`：默认 `False`；`True` 时使用宽松的 `[C,H,S,N,halogen]` 集合。

**返回类型**

- `float`：疏水部分的表面积。

**异常**

- 若 `rdFreeSASA` 或 `Chem.MolFromSmarts` 失败，会直接抛出原始异常。

**注意事项**

- ⚠️ `query` 取 `hydrophobic_query.GetAtoms()[0]`，仅用于指定原子类型；不能传入空 `Smarts`。

**示例**

```python
hydrophobic_area = MolCaculator.calc_hydrophobic_surface_area(mol, is_high_precise=True)
```

---

### `MolCaculator.calc_diameter(mol: Chem.Mol) -> float`

**功能说明**

- 以 3D 原子坐标为基础，利用 `scipy.spatial.distance.pdist` 获取所有原子对中最远距离。

**参数**

- `mol`：必须有至少一个 3D conformer。

**返回类型**

- `float`：最大直径（Å）。

**异常**

- `ValueError`：若 `mol.GetNumConformers() == 0`。

**注意事项**

- 当原子数 < 2 时直接返回 `0.0`；多原子建议先调用 `Chem.AddHs` 以考虑氢原子对距离的影响。

**示例**

```python
diam = MolCaculator.calc_diameter(mol)
```

---

### `MolCaculator.calc_sssr(mol: Chem.Mol) -> Dict[str, Any]`

**功能说明**

- 结构化地返回 SSSR、RingInfo 的原子/键环、原子/键环计数及 `SymmSSSR` 算法名。

**参数**

- `mol`：任意 RDKit `Mol`；可以不包含 3D。

**返回类型**

- `Dict`：包含 `"num_rings"`、`"atom_rings"`、`"by_size"`、`"atom_ring_count"` 等字段。

**异常**

- `ValueError`：若 `mol` 为 `None` 或 `Chem.GetSymmSSSR` 解析失败。

**注意事项**

- 💡 `ri_atom_rings`/`ri_bond_rings` 可能与 SSSR 不完全一致；`atom_in_ring/bond_in_ring` 来自计数是否 > 0。

**示例**

```python
sssr_info = MolCaculator.calc_sssr(mol)
```

---

### `MolCaculator.calc_logp(mol: Chem.Mol) -> float`

**功能说明**

- 基于 Crippen 模型计算 LogP。

**参数**

- `mol`：任意 RDKit 分子。

**返回类型**

- `float`：化合物的 LogP。

**异常**

- 纯函数，无显式异常。

**注意事项**

- 💡 不会自动处理显式氢；若需要氢参照 `Chem.AddHs`。

**示例**

```python
logp = MolCaculator.calc_logp(mol)
```

---

### `MolCaculator.calc_morgan_fringerprint(mol: Chem.Mol, countSimulation: bool = False, bitSize: int = 2048) -> np.array`

**功能说明**

- 生成带计数选项的 Morgan 指纹，最终转换为 NumPy 数组（浮点）。

**参数**

- `mol`：RDKit 分子。
- `countSimulation`：开启后指纹包含计数信息，否则为二进制指纹。
- `bitSize`：指纹维度。

**返回类型**

- `np.array`：单行数组，长度与 `bitSize` 匹配。

**异常**

- `DataStructs.ConvertToNumpyArray` 上游异常。

**注意事项**

- ⚠️ 返回值是 `copy.deepcopy(arr)`，确保不会与 RDKit 内部缓冲共享；`countSimulation=True` 适合需要计数特征的模型。

**示例**

```python
fingerprint = MolCaculator.calc_morgan_fringerprint(mol, countSimulation=True, bitSize=4096)
```

---

### `MolCaculator.calc_maccs_fingerprint(mol: Chem.Mol) -> np.array`

**功能说明**

- 调用 `MACCSkeys.GenMACCSKeys` 生成 166 位指纹并以整数数组返回。

**参数**

- `mol`：所需分子对象。

**返回类型**

- `np.array`：形状 `(1,)` 的整数数组，内部值为 0/1。

**异常**

- RDKit 计算指纹异常会贯穿到调用方。

**注意事项**

- 💡 默认 `dtype=np.int32`，如果在训练管线需 float 可自行转。

**示例**

```python
maccs = MolCaculator.calc_maccs_fingerprint(mol)
```

---

### `MolCaculator.calc_crippen_contribs(mol: Chem.Mol) -> tuple`

**功能说明**

- 通过 `rdMolDescriptors._CalcCrippenContribs` 返回每个原子的 LogP 与 MolRefr 贡献统计。

**参数**

- `mol`：建议事先调用 `mol.UpdatePropertyCache()` 和 `AllChem.ComputeGasteigerCharges(mol)`（已在方法内执行）。

**返回类型**

- `tuple`：通常包含两个列表，分别为平均贡献与标准差。

**异常**

- 若 RDKit 计算失败将抛出原始异常。

**注意事项**

- ⚠️ `_CalcCrippenContribs` 属于内部 API，若 future RDKit 版本调整需留意。

**示例**

```python
crippen = MolCaculator.calc_crippen_contribs(mol)
```

---

### `MolCaculator.get_ring_atoms(mol, use_ringinfo: bool = True) -> List[List[int]]`

**功能说明**

- 返回所有环的原子索引列表，默认用 `ringinfo`；可切换为 `GetSymmSSSR` 做兼容。

**参数**

- `mol`：RDKit 分子。
- `use_ringinfo`：`True` 时访问 `mol.GetRingInfo()`，`False` 时退回 `GetSymmSSSR`。

**返回类型**

- `List[List[int]]`：每个环内原子索引（0-based）。

**异常**

- 无。

**注意事项**

- 💡 注释中提到可转为 1-based，当前实现保留 0-based 以和其余函数一致。

**示例**

```python
rings = MolCaculator.get_ring_atoms(mol)
```

---

### `MolCaculator.get_hydrophobic_clusters(mol: Chem.Mol, distance_threshold: float = 2.0, extended: bool = True) -> List[List[int]]`

**功能说明**

- 在 3D 空间中以距离阈值为依据，把满足邻居规则的碳原子聚成簇；返回 1-based 索引。

**参数**

- `mol`：要求有 3D 构象；内部会 embed+optimize（若尚未嵌入）。
- `distance_threshold`：聚类阈值。
- `extended`：决定邻居元素判断。

**返回类型**

- `List[List[int]]`：每个子列表代表一个疏水簇。

**异常**

- `math`、RDKit embed 失效会抛出；当前调用中没有额外捕获。

**注意事项**

- ⚠️ 即使 `mol` 没有 3D 构象也会尝试 `AllChem.EmbedMolecule`，但 embed 失败不会抛出而继续处理。

**示例**

```python
clusters = MolCaculator.get_hydrophobic_clusters(mol, distance_threshold=1.5)
```

---

### `MolCaculator.approximate_oe_hydrophobe(...)`

**功能说明**

- 近似 OEShape 的疏水簇识别，引入部分电荷、芳香性、最小簇大小等多重过滤；返回 `(count, [atom_ids])` 元组列表。

**参数**

- `mol`：可选自动 embed；若 `embed_if_needed` 且无 conformer，函数会尝试 `ETKDGv3` + `MMFF`。
- 其余参数控制阈值（详见代码），如 `partial_charge_cutoff`、`min_cluster_size`、`extended_filter`。

**返回类型**

- `List[Tuple[int, List[int]]]`：按 `(簇大小, [1-based 原子索引])` 排序。

**异常**

- Embed/优化失败时仅打印，不抛；其余异常会冒泡。

**注意事项**

- ⚠️ `partial_charge_cutoff` 需配合 Gasteiger 计算，会修改 `mol` 的 `_GasteigerCharge` 属性；`embed_if_needed=False` 时请确保己有 3D。

**示例**

```python
hydrophobes = MolCaculator.approximate_oe_hydrophobe(mol, distance_threshold=1.2)
```

---

### `MolCaculator.get_anions_cations(mol) -> Tuple[List[int], List[int]]`

**功能说明**

- 按 `atom.GetFormalCharge()` 分类，分别返回阴离子与阳离子原子索引（0-based）。

**参数**

- `mol`：任意 `Chem.Mol`。

**返回类型**

- `Tuple[List[int], List[int]]`。

**异常**

- 无。

**注意事项**

- 💡 对于多价电荷原子也会被正确捕捉。

**示例**

```python
anions, cations = MolCaculator.get_anions_cations(mol)
```

---

### `MolCaculator.get_hbond_acceptors(mol)`

**功能说明**

- 通过 SMARTS 匹配常见受体，返回 `(atom_idx, ...)` 的 tuple 列表。

**参数**

- `mol`：需要的 RDKit 分子。

**返回类型**

- `Tuple[Tuple[int, ...], ...]`（`GetSubstructMatches` 输出）。

**异常**

- 若 SMARTS 语法无效，RDKit 会抛出 `Chem.Atom.GetSubstructMatches` 异常。

**注意事项**

- ⚠️ SMARTS 模式包含 `O`/`n` 等，若需要更具体识别可自行替代。

**示例**

```python
acceptors = MolCaculator.get_hbond_acceptors(mol)
```

---

### `MolCaculator.get_hbond_donors(mol)`

**功能说明**

- 同样通过 SMARTS 获取氢键供体元组。

**参数**

- `mol`：RDKit 分子。

**返回类型**

- `Tuple[Tuple[int, ...], ...]`

**异常**

- 同上。

**注意事项**

- 💡 可配合 `get_hbond_acceptors` 生成受体/供体的标注。

**示例**

```python
donors = MolCaculator.get_hbond_donors(mol)
```

---

### `MolCaculator.calc_pharmacophore_features(mol: Chem.Mol, extended_hydrophobe: bool = True, use_ringinfo: bool = True)`

**功能说明**

- 聚合上述各模块，提取氢键受体/供体、离子、环、疏水簇，返回特征字典（存在性/计数/原子列表）。

**参数**

- `mol`：无 conformer 时会尝试 embed + UFF。
- `extended_hydrophobe`：控制疏水聚类使用的过滤严格程度。
- `use_ringinfo`：决定环原子提取路径。

**返回类型**

- `Tuple[Dict[str, int], Dict[str, int], Dict[str, Any]]`。

**异常**

- 所有子方法异常会冒泡，但 embed/optimize 的失败被吞掉。

**注意事项**

- ⚠️ 返回 `feature_atom_list` 中的疏水簇对应 `approximate_oe_hydrophobe`，即 1-based。

**示例**

```python
features, counts, atoms = MolCaculator.calc_pharmacophore_features(mol)
```

---

### `MolCaculator.calc_intra_pharmacophore_distance(mol: Chem.Mol, atom_list: dict, conf_id: int = 0)`

**功能说明**

- 依照给定的 atom_list（1-based）计算每种药效团内部原子对之间的距离矩阵。

**参数**

- `mol`：必须有指定 `conf_id` 的 3D conformer。
- `atom_list`：键为药效团类型，值为原子索引列表。
- `conf_id`：默认 0。

**返回类型**

- `Dict[str, List[List[float]]]`。

**异常**

- `ValueError`：若 conformer 缺乏 3D 坐标。

**注意事项**

- ⚠️ 假定 `atom_list` 中的索引为 1-based，会在内部减 1。

**示例**

```python
distances = MolCaculator.calc_intra_pharmacophore_distance(mol, atoms_dict)
```

---

### `MolCaculator.calc_inter_pharmacophore_distance(mol, atom_list, type1, type2, conf_id: int = 0) -> List[List[float]]`

**功能说明**

- 计算两种药效团类型之间的距离矩阵，返回空列表表示任一类型无对应原子。

**参数**

- `mol`：含目标 conformer。
- `atom_list`：同上。
- `type1`, `type2`：例如 `'rings'`、`'acceptor'`。
- `conf_id`：目标 conformer。

**返回类型**

- `List[List[float]]`。

**异常**

- `ValueError`：若 conformer 缺 3D。

**注意事项**

- 💡 可用于计算离群药效团间距以供可视化或聚类。

**示例**

```python
inter = MolCaculator.calc_inter_pharmacophore_distance(mol, atoms_dict, "donor", "acceptor")
```
