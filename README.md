# AI Lingues Biotech Python Library

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/status-active-brightgreen.svg)](https://github.com/ailingues/pybiotech)


pybiotech is a biology-focused toolkit from AI Lingues that builds on RDKit, PubChem, and UniProt to deliver molecule optimization, fingerprinting, structured loaders, and request modeling for data-driven biomedical workflows.

---

## Key Capabilities

- Molecular engineering with ETKDG embedding, MMFF/UFF optimization, surface and hydrophobic metrics, ring counts, fingerprints, and pharmacophore distances.
- Binary serialization that wraps rdkit.Chem.Mol objects in a CRC32-protected container suitable for IPC, sockets, or shared memory.
- Loader stack covering file, directory, and text SDF parsing plus UniProt XML streaming with error tolerance and namespace handling.
- PubChem request builder that unifies Input, Operation, Output, and Query segments, enforces validation rules, and can reverse-parse existing URLs.
- Schema modeling using Pydantic and xsdata so UniProt query fields and PubChem domain pieces stay type-safe.

---

## 📦 Installation

### Install from PyPI

```bash
pip install pybiotech
```

### Requirements

- Python 3.10 or newer (3.11 recommended)
- RDKit for molecule processing
- Additional dependencies declared in pyproject.toml / requirements.txt (xsdata, lxml, pydantic, scipy, numpy, etc.)

---

## 📚 Quick Start

### 1. Embed, optimize, and serialize a molecule

```python
from rdkit import Chem
from pybiotech.core.molecule.optimizer import Optimizer

mol = Chem.MolFromSmiles("c1ccccc1O")
mol3d, ok, meta = Optimizer.embed_and_optimize(mol, random_seed=42)
blob = Optimizer.to_serialize([mol3d, None], include_props=True, with_checksum=True)
print("converged:", ok, "method:", meta["method"])
```

### 2. Read SDF records from a directory

```python
from pybiotech.loaders.sdf_loader import SDFLoader

mols, expected, actual = SDFLoader.readDataFromDir("data/sdf", recursive=True)
print(f"expected {expected}, loaded {actual} valid molecules")
```

### 3. Build a PubChem REST URL

```python
from pybiotech.classes.nih.https.pubchem.input_spec import InputSpec
from pybiotech.classes.nih.https.pubchem.opera_spec import OperationSpec
from pybiotech.classes.nih.https.pubchem.output_spec import OutputSpec
from pybiotech.classes.nih.https.pubchem.query_options import QueryOptions
from pybiotech.classes.nih.https.pubchem.request import PubChemRequest

input_spec = InputSpec(domain="compound", namespace="cid", identifiers="2244")
operation_spec = OperationSpec(operation="property", tags=["MolecularFormula", "InChIKey"], domain="compound")
output_spec = OutputSpec(output_format="JSON")
query_options = QueryOptions(record_type="3d", image_size="large")

req = PubChemRequest(
    input_spec=input_spec,
    operation_spec=operation_spec,
    output_spec=output_spec,
    query_options=query_options,
)
print(req.build_url())
```

---

## Module overview

### Core (molecule)

| Module | Purpose | Docs |
| --- | --- | --- |
| `core.molecule.calculator` | Geometry, surface, fingerprint, and pharmacophore helpers for Chem.Mol. | [docs/core/molecule/calculator.md](docs/core/molecule/calculator.md) |
| `core.molecule.optimizer` | ETKDG embedding plus MMFF/UFF optimization, serialization container, and checksum-aware deserialization. | [docs/core/molecule/optimizer.md](docs/core/molecule/optimizer.md) |

### Loaders

| Module | Purpose | Docs |
| --- | --- | --- |
| `loaders.sdf_loader` | Stream file/dir/text SDF content with index filtering, error tolerance, and warning control. | [docs/pybiotech/loaders/sdf_loader.md](docs/pybiotech/loaders/sdf_loader.md) |
| `loaders.uniprot_loader` | Auto namespace detection and entry-by-entry UniProt XML parsing. | [docs/pybiotech/loaders/uniprot_loader.md](docs/pybiotech/loaders/uniprot_loader.md) |
| `loaders.uniprot_query_field_loader` | Load QueryField definitions from JSON/text and build validated term:value fragments. | [docs/pybiotech/loaders/uniprot_query_field_loader.md](docs/pybiotech/loaders/uniprot_query_field_loader.md) |

### Classes (models and request helpers)

| Module | Purpose | Docs |
| --- | --- | --- |
| `classes.uniprot.https.uniprot.org.uniprot_query_field` | Pydantic QueryField model with recursive siblings/items. | [docs/classes/uniprot/https/uniprot/org/uniprot_query_field.md](docs/classes/uniprot/https/uniprot/org/uniprot_query_field.md) |
| `classes.uniprot.https.uniprot.org.uniprot` | xsdata-generated dataclasses and enums for UniProt Entry, Protein, Sequence, etc. | [docs/classes/uniprot/https/uniprot/org/uniprot.md](docs/classes/uniprot/https/uniprot/org/uniprot.md) |
| `classes.nih.https.pubchem.*` | Input, Operation, Output, Query, Request, and Identifier/domain modeling for PubChem PUG REST. | [docs/classes/nih/https/pubchem](docs/classes/nih/https/pubchem) |

---

## 📖 Documentation

Detailed API references are under docs/:

- docs/core/molecule/ - optimizer and calculator references
- docs/pybiotech/loaders/ - SDF/UniProt loader documentation
- docs/classes/nih/https/pubchem/ - PubChem URL and query modeling
- docs/classes/uniprot/https/uniprot/org/ - UniProt schema reference

---

## 🏗️ Project layout

- pybiotech/
  - pybiotech/
    - core/
      - molecule/
    - loaders/
    - classes/
  - docs/
  - scripts/
  - README.md
  - README_Example.md
  - pyproject.toml
  - requirements.txt

---

## 📝 License & Usage Terms

### Package Usage License

The PyBiotech package (the binary package installed via pip) is released under the MIT License:

- ✅ Free Use – May be used freely in personal and enterprise projects
- ✅ Commercial Use – Allowed in commercial products and services
- ✅ Free Distribution – You may freely distribute and redistribute the package
- ✅ No Usage Restrictions – No fees or additional authorization required

### Source Code Protection Terms

Important Notice: The source code of this project is private property and protected by intellectual property law:

- ❌ Source Code Not Public – The source code is not publicly available
- ❌ No Reverse Engineering – Decompilation, reverse engineering, or disassembly of the package is strictly prohibited
- ❌ No Source Distribution – You may not obtain, copy, or distribute the source code in any form
- ❌ No Modified Redistribution – You may not modify the package and redistribute it
In short: You are free to use our package (including for commercial purposes), but please respect our source code intellectual property rights.

## 📧 Contact Us

Website: https://www.ailingues.com

Email: support@ailingues.com

Technical Support: For any questions or suggestions, please contact us via email.

---

**Made with ❤️ by AI Lingues Team** 
=======
- 支持NIH PubChem公开数据库数据查询访问
- 支持SDF格式文件高速读取
- 支持分子化合物构象力场优化及特征计算

<hr>

**最新版本**

version **0.2.6**

**主要依赖**

- python >= 3.11
- pycorelibs >= 0.2.6
- rdkit == 2024.9.6

**CopyRight**

    AI Lingues Team

**email**

    support@ailingues.com

<hr>

# core 模块

<hr>

## molecule 模块

### calculator 模块

<hr>

#### MolCaculator 分子计算器类

计算分子化合物相关的参数、特征等，函数清单如下：

| 函数                              | 描述                                                                   |
| --------------------------------- | ---------------------------------------------------------------------- |
| approximate_oe_hydrophobe         | 计算"疏水簇"                                                           |
| calc_surface_area                 | 计算分子化合物的表面积                                                 |
| calc_hydrophobic_surface_area     | 计算分子疏水表面积                                                     |
| calc_diameter                     | 计算单个分子的最大直径                                                 |
| calc_sssr                         | 计算分子的环信息(SSSR 与 RingInfo 快照）                               |
| calc_logp                         | 计算分子化合物的脂溶性(LogP)                                           |
| calc_morgan_fringerprint          | 计算分子化合物的Morgan指纹                                             |
| calc_maccs_fingerprint            | 计算分子化合物的MACCS指纹                                              |
| calc_crippen_contribs             | 计算分子的Crippen规范化贡献,包括每个原子的 logP 贡献和摩尔折射率贡献   |
| calc_pharmacophore_features       | 计算药效团特征                                                         |
| calc_intra_pharmacophore_distance | 计算每种药效团类型内部所有原子的欧几里得距离矩阵                       |
| calc_inter_pharmacophore_distance | 计算两个不同药效团类型之间的欧几里得距离矩阵                           |
| get_ring_atoms                    | 提取所有环中原子的索引 (1-based)                                       |
| get_hydrophobic_clusters          | 识别疏水原子并根据 3D 距离聚合成多个“簇”，以对应 SDF 中多行 hydrophobe |
| get_anions_cations                | 获取阴离子原子列表、阳离子原子列表 (1-based)                           |
| get_hbond_acceptors               | 识别氢键受体原子 (0-based 索引)                                        |
| get_hbond_donors                  | 识别氢键供体原子 (0-based 索引)                                        |

<hr>

##### approximate_oe_hydrophobe

计算"疏水簇"。  
在纯 RDKit 环境下近似地模仿 OEShape 的疏水原子聚合逻辑, 返回多个"疏水簇"。
每个簇可视为一行 hydrophobe, 类似:  
3  38 41 42 hydrophobe  
4  17 19 20 21 hydrophobe

等等.

###### 参数说明

- mol : RDKit Mol 对象

    如果已带有合理的 3D conformer，可不再 embed。如果没3D,且embed_if_needed=True则会自动Embed+Optimize。

- distance_threshold : float

    两个候选原子若距离 < 此阈值则视为同一疏水簇. OEShape常用 ~1.0 or 1.5Å.

- exclude_aromatic : bool

    是否排除所有芳香碳(例如苯环C). OEShape 里通常某些芳香C也可能算疏水, 这里可选。

- partial_charge_cutoff : float or None

    若不为 None, 则使用 Gasteiger 部分电荷, 排除绝对值>=该阈值的碳, 以排除极性碳.  
    例如 0.2 -> |q|≥0.2 的碳视为不疏水.

- min_cluster_size : int

    最小簇大小, 若某个簇只有 < min_cluster_size 个原子, 则可视为噪声/舍弃 (或可保留).

- extended_filter : bool

    若 True, 用 "不含 O,N,S,P,卤素" 规则剔除碳; 否则只要邻居里无 O,N 即保留。

- embed_if_needed : bool

    若 mol 无 3D 构象, 是否调用 ETKDG embed.

- max_attempts : int

    embed 出错时的尝试次数.

###### 返回

- hydrophobe_lines : list of tuples  [(atom_count, [a1,a2,...]), ...]

    其中 a1,a2,... 是1-based原子索引, 已排序, 代表同一个疏水簇.
    可以把它转成 SDF-like字符串:
        f"{atom_count} {' '.join(map(str, atom_ids))} hydrophobe"
    也可以直接当数据结构用.

###### 注意

1. 不保证与 OEShape 结果完全一致, 只是"力场 + 距离聚类 + 不同过滤"的思路.
2. 若 embed 失败/坐标不合理, 或分子过大, 可能结果仍不理想.
3. 可多次调参 distance_threshold, partial_charge_cutoff 等, 观察对结果的影响.

##### calc_surface_area

计算分子化合物的表面积。

在调用之前,必须做如下处理：

 1. 添加显式氢原子
    ```python
        mol.UpdatePropertyCache(strict=False)  
        mol_with_H = Chem.AddHs(mol)
    ```

 2. 生成 3D 坐标
    ```python
        AllChem.EmbedMolecule(mol_with_H)  
        AllChem.MMFFOptimizeMolecule(mol_with_H)
    ```

###### Args

- mol (Chem.Mol): 分子化合物对象

###### Raises

- Exception: _description_

###### Returns

- float: 分子总表面积

<hr>

##### calc_hydrophobic_surface_area

计算分子疏水表面积。

在调用之前,必须做如下处理：

1. 添加显式氢原子
    ```python
        mol.UpdatePropertyCache(strict=False)  
        mol_with_H = Chem.AddHs(mol)
    ```

2. 生成 3D 坐标
    ```python
        AllChem.EmbedMolecule(mol_with_H)  
        AllChem.MMFFOptimizeMolecule(mol_with_H)
    ```

###### Args

- mol (Chem.Mol): 分子化合物对象
- is_high_precise (bool): 是否更精确,缺省为False

###### Raises

- e: _description_

###### Returns

- float: 分子疏水表面积

<hr>

##### calc_diameter

计算单个分子的最大直径。

###### Args

- mol (Chem.Mol): RDKit分子对象,假设已经有3D坐标。

###### Returns

- float: 分子的最大直径（Angstrom）。

<hr>

##### calc_sssr

计算分子的环信息（SSSR 与 RingInfo 快照）。

###### 功能

- 调用 GetSymmSSSR(mol) 触发并获取“对称 SSSR”环集（以原子索引表示）。
- 读取 RingInfo（mol.GetRingInfo()），给出每个原子/键的“环计数”等信息，以及按原子/按键的环列表。
- 返回结构化结果，便于后续分析与统计（不再 print）。

###### 参数

- mol : rdkit.Chem.rdchem.Mol

    RDKit 分子对象。函数内部不会修改该对象（仅读取）。

###### 返回

- Dict[str, Any]
    {

    "num_rings": int,                     # SSSR 环的数量

    "atom_rings": List[List[int]],        # SSSR：每个环对应的原子索引列表

    "by_size": Dict[int, List[List[int]]],# 按环尺寸分组的 SSSR

    "ri_atom_rings": List[List[int]],     # RingInfo.AtomRings()（不一定等同于 SSSR）

    "ri_bond_rings": List[List[int]],     # RingInfo.BondRings()

    "atom_ring_count": List[int],         # 每个原子属于多少个环

    "bond_ring_count": List[int],         # 每根键属于多少个环

    "atom_in_ring": List[bool],           # 原子是否在任何环中（派生自计数>0）

    "bond_in_ring": List[bool],           # 键是否在任何环中（派生自计数>0）

    "algorithm": str,                     # 'SymmSSSR'

    }

###### 说明

- GetSymmSSSR() 会确保环感知已进行，并把信息缓存到 RingInfo。
- 返回中的 `atom_rings` 是 SSSR（对称最小环集）；`ri_atom_rings/ri_bond_rings` 来自 RingInfo，可能包含与 SSSR 不完全一致的环枚举（实现层面差异）。
- 原子/键的“是否在环中”通过计数 > 0 派生，效率高且直观。

<hr>

##### calc_logp

计算分子化合物的脂溶性(LogP)。

###### Args

- mol (Chem.Mol): 分子化合物Mol实例对象,支持Chem.Mol子类实例

###### Returns

- float: 脂溶性(LogP)值

<hr>

##### calc_morgan_fringerprint

计算分子化合物的Morgan指纹。

###### 说明

Morgan指纹是RDKit中一种常用的分子指纹类型,可以用于描述分子的结构和相似性。

它基于分子的拓扑结构和半径参数生成,具有以下特点：

1. 生成的指纹是一个固定长度的二进制向量,每个位表示一个子结构的存在或缺失。
2. 指纹的长度和半径参数可以根据需要进行调整,以平衡指纹的信息量和计算效率。
3. 可以使用不同的哈希函数来生成指纹,以增加指纹的多样性和鲁棒性。

GetMorganGenerator签名: 参考doc/specifications/interface/GetMorganGenerator.md

###### 注意事项

关于countSimulation参数

1. Morgan 指纹默认行为:

    - 默认情况下（countSimulation=False）：

        Morgan 指纹是一个位向量,值为 0 或 1,表示某个化学环境是否存在。

    - 启用计数模拟（countSimulation=True）：

        Morgan 指纹包含整数值,表示某个化学环境出现的次数。

2. 在分类问题中：

    2.1 如果化学环境的 存在与否 是关键,则 0 和 1 的位向量形式通常足够。

      - 适用场景：

        a.化学环境的存在与否足够描述目标性质。

        b.任务是分类问题（例如,是否具有毒性、是否活跃）。

        c.数据稀疏,或子结构的出现次数分布较均匀。

    2.2 如果化学环境的 出现频率 是分类的潜在决定因素,则保留计数信息可能更有帮助。

      - 适用场景：

        a.化学环境的出现频率对分类任务有重要影响。

        b.任务需要描述分子中功能性团的强度（如高毒性分子）。

        c.需要捕捉数量信息的额外价值。

        d.其他解决预测问题或者分析场景

###### 引用[fingerprint](https://github.com/daiyizheng/DL/blob/master/07-rdkit/08-rdkit%E5%8C%96%E5%AD%A6%E6%8C%87%E7%BA%B9.ipynb)

###### Args

- mol (Chem.Mol): 分子化合物Mol实例对象,支持Chem.Mol子类实例
- countSimulation (bool): 是否开启计数,缺省为False(此参数详细参考注意事项部分)
- bitSize (int): 位向量长度,缺省2048

###### Returns

- np.array: Morgan指纹数据数组

<hr>

##### calc_maccs_fingerprint

计算分子化合物的MACCS指纹。

###### 方法

使用rdkit.Chem.MACCSkeys.GenMACCSKeys 函数来计算分子

###### 说明

MACCS (Molecular ACCess System) 分子指纹是一种用于表示分子结构信息的二进制指纹。  
MACCS分子指纹是基于分子中是否含有特定的亚结构来定义的,共包含166个不同的分子特征。  
每个特征都对应于一个特定的化学子结构,例如,一个羟基、一个苯环或一个氮原子等。  
如果分子中存在这个特征,则该特征对应的二进制位上的值为1,否则为0。  
MACCS分子指纹的长度为166位,它可以用于分子相似性比较、分子分类、分子聚类、分子筛选等许多领域中的化学信息学研究。

###### 注意事项

无

###### 引用 [fingerprint](https://github.com/daiyizheng/DL/blob/master/07-rdkit/08-rdkit%E5%8C%96%E5%AD%A6%E6%8C%87%E7%BA%B9.ipynb)

###### Args

- mol (Chem.Mol): 分子化合物Mol实例对象,支持Chem.Mol子类实例

###### Returns

- np.array: MACCS指纹数据数组

<hr>

##### calc_crippen_contribs

计算分子的Crippen规范化贡献,包括每个原子的 logP 贡献和摩尔折射率贡献。

###### 方法

无

###### 说明

基于分子的原子电荷和分子的几何形状计算的,可以用于描述分子的溶解度、生物利用度和其他性质。  
这个函数通常与RDKit分子对象一起使用。

###### 注意事项

  1. 传入的Chem.Mol对象应先调用UpdatePropertyCache方法处理
  2. Crippen规范化贡献虽然是按整个分子化合物计算,但计算结果应按索引位置将贡献值分配到对应原子作为原子特征的一部分

###### Args

- mol (Chem.Mol): 分子化合物Mol实例对象,支持Chem.Mol子类实例

###### Returns

- tuple: 包括每个原子的 logP 贡献和摩尔折射率贡献
        元组,其中包含两个长度为分子中原子数的列表。
        第一个列表包含每个原子的Crippen贡献的平均值,
        第二个列表包含每个原子的Crippen贡献的标准差。

<hr>

##### calc_pharmacophore_features

计算药效团特征。

###### 说明

  1. 检测是否含 3D conformer, 若无则做简单的Embed + Optimize(可选).
  2. 计算氢键受体/供体, 阴阳离子, 环原子, 疏水原子等.
  3. 返回 (features, features_count, atom_list).
  其中:
      - features = {"rings":0/1,...}
      - features_count = {"rings":N,...}
      - atom_list = {"rings":[...],...} (1-based or list of lists)

<hr>

##### calc_intra_pharmacophore_distance

计算每种药效团类型内部所有原子的欧几里得距离矩阵。

###### 参数

- mol: RDKit Mol 对象,需包含 3D 坐标 (Conformer)。
- atom_list: dict,每种药效团类型对应的原子编号列表,如 {'rings': [1,2,3], 'anion': [4,5], ...}
- conf_id: int,可选,指定使用哪个 conformer 计算距离。

###### 返回

- intra_distances: dict

    key为药效团类型,value为对应的距离矩阵 (二维list),  
    如: {'rings': [[0.0, 1.2, ...], [...], ...], 'anion': [...], ...}

<hr>

##### calc_inter_pharmacophore_distance

计算两个不同药效团类型之间的欧几里得距离矩阵。

###### 参数

- mol: RDKit Mol 对象,需包含 3D 坐标 (Conformer)。
- atom_list: dict,每种药效团类型对应的原子编号列表,例如:

    {'rings': [1,2,3], 'anion': [4,5], 'cation': [], ...}

- type1: str,第一个药效团类型 (如 'rings', 'anion', 'cation', 'acceptor', 'donor', 'hydrophobe')
- type2: str,第二个药效团类型
- conf_id: int,可选,指定使用哪个 conformer 计算距离。

###### 返回

- inter_distance_matrix:

    二维 list, 形状为 (len(type1原子), len(type2原子))

<hr>

##### get_ring_atoms

提取所有环中原子的索引 (1-based)。

###### 参数

- mol: RDKit Mol 对象
- use_ringinfo: bool

    如果为 True, 使用 ringinfo 来识别环原子;  
    如果为 False, 使用 GetSymmSSSR.

###### 返回

- List[List[int]] : 每个环是一个列表, 里面存环内的原子(1-based).

    例如: [[1,2,3,4,5,6],[8,9,10]].

<hr>

##### get_hydrophobic_clusters

识别疏水原子并根据 3D 距离聚合成多个“簇”，以对应 SDF 中多行 hydrophobe。

返回: List[List[int]], 每个子列表是一群(簇)疏水原子的 1-based 索引。

###### 参数

- mol : RDKit Mol 对象 (需有3D构象,若无需先 Embed + 优化)
- distance_threshold : float

    任意两个候选疏水原子的3D距离若 < 该值，就视为同一簇。  
    默认1.0Å，也可尝试1.5/2.0等。
- extended : bool

    True: 不含 O,N,S,P,卤素(F,Cl,Br,I)的碳视为疏水  
    False: 仅要求邻居里无 O,N

###### 返回

- clusters_1based : List[List[int]]

    例如 [[10,12,14],[18,22,23]]，表示两簇疏水原子(1-based)。  
    若没有疏水原子，返回空列表 []。

<hr>

##### get_anions_cations

获取阴离子原子列表、阳离子原子列表 (1-based)。

当前基于 formal charge 判定:  
- atom.GetFormalCharge() <0 => anion
- atom.GetFormalCharge() >0 => cation

对多价电荷, 同样识别到同一组, 如 +2 => cation.  
若需部分电荷, 需额外力场/量化计算.

<hr>

##### get_hbond_acceptors

识别氢键受体原子 (0-based 索引)。

返回 SubstructMatch 的 tuple list, 每个元素是 (atom_idx, ...).  
如果只需原子 idx, 可自行提取 match[0].  
这里使用稍微更全的 SMARTS 例子, 包含芳环N, 羰基O等.

<hr>

##### get_hbond_donors

识别氢键供体原子 (0-based 索引)。

<hr>

### optimizer 模块

分子构象优化，序列化/反序列化

#### Optimizer 构象优化类

<hr>

##### embed_and_optimize (静态类方法)

使用RDKit生成初始3D构象并优化几何结构。优先使用MMFF94力场优化，失败则回退到UFF。

###### Args

- mol (Mol): RDKit 分子对象（可未消毒）。函数内部会复制一份工作副本，不会修改来参。
- max_embed_attempts (int, optional): 3D 构象嵌入（ETKDG）的最大尝试次数。

    数值越大，困难分子的成功率越高，但时间也越长。  
    `建议`：一般 200–1000；含大环/复杂稠环可适当提高。. Defaults to 1000.
- random_seed (int, optional): 随机种子。固定值可复现结果；

    设置为 -1 表示完全随机（非确定性）。  
    `建议`：科研/调试阶段建议固定；生产批量可使用 -1 提高多样性. Defaults to 0xC0FFEE.
- use_small_ring_torsions (bool, optional): ETKDG 的小环扭转参数。

    开启通常更符合小环（如 3–5 元环）经验构象，提升嵌入质量。 to True.
- use_macrocycle_torsions (bool, optional): . ETKDG 的大环扭转处理。

    对大环/多环体系开启有利于找到更合理的初始构象。 to True.
- prune_rms_thresh (float, optional):

    构象剔除的 RMSD 阈值（重复构象的去冗策略）。  
    即便只嵌 1 个构象，这个阈值也会影响“寻找与已有构象足够不同”的重试逻辑。  
    值越大，越容易把相似构象视为“重复”而继续尝试。  
    `建议`：0.1–0.5 Å 之间较常用。. Defaults to 0.1.
- max_ff_iters (int, optional): 力场最小化的最大迭代次数。

    数值越大，越有机会“收敛”；但时间也更长。

    `建议`：200–1000。若经常“不收敛”，可先增大再考虑结构预处理。. Defaults to 500.

###### **Raises:**

- ValueError: RDKit 分子对象为空

###### **Returns:**

Tuple[ Mol, bool, Dict[str, Any]]:

- optimized_mol : rdkit.Chem.Mol 已添加显式氢的分子对象，包含单一 3D 构象。

    若过程中失败，也会返回当前工作副本以便诊断。
- ok : bool
    是否达到力场收敛条件（True=收敛；False=未收敛/失败）。
- meta : Dict[str, Any]
    诊断信息字典，常见键如下（按情况部分缺省）：
  - stage : str
      当前执行阶段："init" | "sanitize" | "embed" | "optimize"。
  - method : str
      实际使用的力场方法："MMFF94" 或 "UFF"。
  - energy : float
      最终力场能量（力场单位，通常可视为 kcal/mol；仅在同一力场内比较具有可比性）。
  - steps : int
      _RDKit Minimize 返回码_（注意：不是实际步数）。0 表示收敛，非 0 表示未收敛。
  - message : str
      提示/警告/错误信息（例如 "ETKDG embedding failed"、"MMFF params unavailable, fallback to UFF"）。

###### 行为与保证

- 不修改传入的 `mol`；在其复制体上操作。
- 执行消毒（Sanitize）与立体化学分配（AssignStereochemistry）。
- 添加显式氢（AddHs）。
- 使用 ETKDGv3 进行 3D 嵌入；清空并仅保留 1 个构象。
- 优先尝试 MMFF94；若分子不支持，则回退 UFF。
- `ok=True` 表示力场最小化返回码为 0（达到收敛条件）；否则为 False。
- 发生常见化学问题（嵌入失败、力场不可用等）时不抛异常，而是 `ok=False` 并在 `meta['message']` 给出原因。
    若 `mol is None` 或输入不可用，可能抛出 `ValueError`。

###### 使用建议

- 需要结果可重复：保持固定 `random_seed`。
- 大环/复杂体系：保持 `use_macrocycle_torsions=True`，适当调大 `max_embed_attempts`。
- 经常未收敛：增大 `max_ff_iters`；或先做电荷/价态/金属配位等预处理。
- 批量处理时，建议记录/持久化 `meta`，便于后期追溯与质量筛选（如优先选用 MMFF94 且收敛的结果）

<hr>

##### **to_serialize (静态类方法)**

将一组 RDKit Mol 对象序列化为“单一二进制容器 blob”（高性能、无文本中间态）。

该二进制容器旨在用于**跨进程/跨应用 IPC 或持久化**，完整保留构象、坐标、手性和（可选）属性。

###### **Args:**

- mols : Iterable[Optional[rdkit.Chem.Mol]]

    分子序列；可包含 None（将写出空占位记录，保持位置对应）。
- include_props : Chem.PropertyPickleOptions, default Chem.GetDefaultPickleProperties()

    是否将分子属性（props）一并打包。推荐 True。
- with_checksum : bool, default True

    是否为每条记录附加 CRC32 校验。生产环境强烈建议开启（默认开启）。

###### **Returns:**

- bytes
    自定义二进制容器（v2）。推荐通过管道/Socket/共享内存/文件在进程或应用间传输。

###### **Raises**

- ValueError 输入序列为空。
- RDKit 相关异常,个别分子损坏等导致二进制写出失败时。

###### **性能与可靠性**

- 性能：相对文本（SDF/MolBlock/JSON）通常更小更快。CRC32 为 C 实现，开销很低（每秒 GB 级）。
- 可靠：长度前缀 + CRC32 抵御截断/半包/损坏；大端编码利于跨语言一致性。
- 兼容性：用于在线协作/IPC 非常稳妥；若用于超长期归档，建议额外保留文本/JSON 以防极端跨大版本情况。

<hr>

##### **to_unserialize (静态类方法)**

从 Optimizer.to_serialize() 产出的**二进制容器**还原出 Mol/None 列表（位置一一对应）。

###### **兼容性**

- v2 容器：b"RDKB\\x02"（推荐；含 RDKit 版本与 flags/CRC）
- v1 容器：b"RDKB\\x01"（向后兼容；仅 header + count + [len+payload]，无版本/flags/CRC）
- 裸 RDKit 单体二进制：若魔数不匹配，尝试作为**单体 Mol** 的 RDKit 二进制读取，成功则返回长度为 1 的列表。

###### **Args:**

- serialized_mols : bytes

    二进制容器 blob。

###### **Returns:**

- List[Optional[rdkit.Chem.Mol]]

    解析得到的分子列表；`None` 表示对应位置为空占位（或损坏记录在启用 CRC 下被拒绝）。

###### **Raises**

- ValueError
  - 入参为空；
  - 容器头不合法，且也不是裸 RDKit 二进制；
  - v2/v1 容器数据结构截断或记录长度异常；
  - v2 容器且 CRC 校验不通过（说明数据损坏/被截断/被篡改）。

###### **说明**

- 使用大端解码（network order）。
- v2 容器会读取并忽略 RDKit 版本字符串（可根据需要记录日志/检查兼容）。
- 默认在 v2 下启用 CRC32 校验（如果写端开启了该标志）。

<hr>

# loaders 模块

加载数据模块

## sdf_loader 模块

### SDFLoader 类

SDF格式数据加载器

从SDF格式数据中读取分子化合物数据,支持从文件、从目录和从文本三种方式读取。

<hr>

#### **closeWarning 关闭警告信息**

此方法是禁用rdkit的警告信息输出

##### Args

无

##### Returns

无

<hr>

#### openWarning 打开警告信息

此方法是恢复rdkit的警告信息输出

##### Args

无

##### Returns

无

<hr>

#### readDataFromFile

从指定sdf文件中读取 molecule 数据

##### Args

- sdfDoc (str):

    sdf 文档名（含路径）
- startIndex (int, optional):

    起始索引（包含）. Defaults to 0.
- endIndex (int, optional):

    结束索引（不包含，默认值 -1 代表读取到文件末尾）. Defaults to -1.
- ignore_error (bool, optional):

    是否忽略错误. Defaults to True.

##### **Raises:**

- FileNotFoundError: 指定文件不存在
- ValueError: 开始索引必须是非负数
- ValueError: 结束索引必须大于开始索引，或指定为-1
- ValueError: 数据段错误

##### Returns

- tuple[list[Mol], int, int]:

    已读出 molecule 列表, 应读数量, 实际读取数量

----

#### readDataFromDir

从指定目录中读取所有的sdf文件,并且读取全部文件的molecule数据

##### Args

- sdfDir (str):

    指定读取的文件目录
- recursive (bool, optional):

    是否递归子目录. Defaults to True.
- ignore_error (bool, optional):

    是否忽略错误. Defaults to True.

##### Raises

- FileNotFoundError: 指定目录不存在
- TypeError: 指定目录并非目录类型

##### Returns

- tuple[list[Mol], int, int]:

    已读出molecule列表,应读数量,实际读取数量

<hr>

#### readDataFromText

从给定文本中读取分子化合物数据

##### Args

- sdfText (str):

    包含分子化合物的sdf格式内容的文本
- ignore_error (bool, optional):

    是否忽略错误. Defaults to True.

##### Returns

- tuple[list[Mol], int, int]:

    已读出molecule列表,应读数量,实际读取数量

<hr>

#### splitDataByMarker

根据指定的标记将SDF文本数据拆分为多个分子数据段

##### Args

- sdfText (str):

    包含至少一个分子化合物的SDF格式内容的文本

- marker (str, optional):

    用于拆分的标记字符串. Defaults to "$$$$\n".

- strict (bool, optional):
    
    是否严格模式. Defaults to True.

##### Raises:
- ValueError: 如果严格模式下没有找到标记

##### Returns

- List[str]:

    拆分后的多个分子数据段列表

<hr>

## nih.pubchem.online 模块

在线实时获取NIH PubChem数据

<hr>

### compound模块

Module for fetching PubChem compound records and conformers over the NIH PubChem REST API,
parsing the returned SDF/JSON payloads, and returning structured ALNPCompound / ALNPConformer
objects.

#### Primary responsibilities

- Request compound SDF blocks for a list of PubChem CIDs.
- Optionally request conformer metadata (ConformerID) and then fetch conformer SDFs.
- Parse SDF content via SDFLoader and construct ALNPCompound and ALNPConformer instances.
- Return a mapping from CID (string) to ALNPCompound instances, with conformer data
    attached when requested.

#### Provided function

- **get_compound**(cid_list: List[str], include_conformer: bool = False) -> Dict[str, ALNPCompound]

- **get_similarity_compound**(
    input: EInputType,
    value: str | int,
    operation: EOperationType,
    output: EOutputType,
    threshold: int = 90,
    max_records: int = 10
) -> List[str]:

----

#### get_compound方法

##### **Parameters**

- **cid_list (List[str])**

    Sequence of PubChem CIDs to fetch. Each CID should be convertible to string; callers
    typically pass strings (e.g. ["2244","3672"]) or integers converted to strings.
- **include_conformer (bool, default False)**

    When False, only the primary compound SDF/metadata is fetched and returned.
    When True, the function also queries the conformer metadata endpoint to obtain
    ConformerIDs for each CID, and then requests SDFs for those conformers and attaches
    ALNPConformer objects under each ALNPCompound.

##### **Return value**

- Dict[str, ALNPCompound]
    A dictionary keyed by the CID string. Each value is an ALNPCompound instance
    populated with:  
    - PUBCHEM_COMPOUND_CID (from SDF properties)
    - ROW (raw SDF block text for the molecule)
    - If include_conformer True:
      - CONFORMER_ID: List[str] of conformer IDs reported by PubChem for that CID
      - CONFORMERS: Dict[str, ALNPConformer] keyed by conformer ID; each conformer
          contains `PUBCHEM_CONFORMER_ID`, `PUBCHEM_COMPOUND_CID` and `ROW` (raw conformer SDF).

##### **Behavior and error handling**

- Uses HTTP endpoints:
  - Compound SDF by CID(s): /rest/pug/compound/cid/{cid_list}/SDF?response_type=display
  - Conformer metadata for CIDs: /rest/pug/compound/cid/{cid_list}/conformers/JSON
  - Conformer SDFs by Conformer ID(s): /rest/pug/conformers/{conformer_id_list}/SDF?response_type=display
- Interprets fetch_url(...) return value as a dict with at least "status_code", "success",
    and "content" keys.
- For 200 responses with "success" False, raises ValueError indicating CID(s) not found.
- For HTTP 404 or 503 at top-level requests, raises ValueError with descriptive messages.
- For conformer fetching, any exceptions raised while parsing conformer SDF content are
    caught and printed; partial results may still be returned for other CIDs.
- Logs important events and warnings:
  - Missing molecules in SDF payloads
  - CID mismatches between JSON metadata and SDF properties
  - Conformer IDs being processed

##### **Notes, constraints and assumptions**

- The SDF parsing relies on SDFLoader.splitDataByMarker and SDFLoader.readDataFromText.
    Those functions are expected to return lists of RDKit-like molecule objects (mol.GetProp(...))
    and counts. The code assumes SDF blocks include `PUBCHEM_COMPOUND_CID` and (for conformers)
    `PUBCHEM_CONFORMER_ID` properties.
- The returned dictionary is pre-populated with keys from the input cid_list (strings)
    and values set to ALNPCompound instances only when parsed successfully; entries may remain None
    if the SDF for a given CID could not be parsed or was absent.
- Network reliability and rate limits are outside this module's control; callers should
    handle transient failures or consider retry/backoff when calling get_compound with large lists.

##### **Example**

- Simple usage:

```python

    cids = ["2244", "3672"]
    compounds = get_compound(cids, include_conformer=False)
    # compounds is a dict mapping "2244" -> ALNPCompound(...)

```

- Fetch compounds with conformers:

```python

    compounds_with_confs = get_compound(["2244"], include_conformer=True)
    conf_ids = compounds_with_confs["2244"].CONFORMER_ID
    conformers_map = compounds_with_confs["2244"].CONFORMERS
```

##### **Types referenced**

- ALNPCompound: container/dataclass representing a PubChem compound record and any attached conformer info.
- ALNPConformer: container/dataclass representing a single conformer (ID, parent CID, raw SDF).

Security and privacy

- Requests are made to the public PubChem REST endpoints; no credentials are required.
- Raw SDF content is stored in returned objects' ROW fields and may contain structural or identifier data;
    treat returned data according to your privacy/security policies.

----

#### get_similarity_compound方法

##### **Summary**

Retrieve similar compounds from PubChem and return detailed compound data.
This function builds a PubChem similarity search URL from the provided
input parameters, performs the HTTP request, parses the returned data
(SDF/JSON/TXT), extracts matching PubChem Compound IDs (CIDs), and then
fetches full compound records (and optionally conformers) for those CIDs.

##### Parameters

- **input (EInputType):**

  The type of the search input (for example SMILES, InChI, CID, etc.).

- **value (str | int):**

    The search value. For CID input this must be an integer or a numeric
        string. Value must not be empty.

- **operation (EOperationType):**

  The PubChem operation type (for example CIDS or RECORD).

- **output (EOutputType):**

    The requested output format from PubChem (SDF, JSON, TXT, ...).

- **threshold (int,optional):**

  Similarity threshold (percent). Default is 90.

- **max_records (int, optional):**

    Maximum number of similar records to request. Default is 10.


##### Return value

- List[str]
    A list by the CID string. Each value is string

##### Raises

- ValueError
  - If `value` is empty.
  - If `input` is EInputType.CID and `value` cannot be parsed as an integer.
  - If `operation` is CIDS and `output` is EOutputType.SDF (disallowed).
  - If a RECORD operation requests TXT output (invalid combination).
  - If the HTTP fetch returns a non-200 status or indicates failure.
  - If the returned content is empty when parsing JSON/SDF/TXT responses.

##### Behavior and error handling

- Uses HTTP endpoints:
    - Compound SDF by CID(s): /rest/pug/compound/cid/{cid_list}/SDF?response_type=display
    - Conformer metadata for CIDs: /rest/pug/compound/cid/{cid_list}/conformers/JSON
    - Conformer SDFs by Conformer ID(s): /rest/pug/conformers/{conformer_id_list}/SDF?response_type=display
- Interprets fetch_url(...) return value as a dict with at least "status_code", "success",
    and "content" keys.
- For 200 responses with "success" False, raises ValueError indicating CID(s) not found.
- For HTTP 404 or 503 at top-level requests, raises ValueError with descriptive messages.
- For conformer fetching, any exceptions raised while parsing conformer SDF content are
    caught and printed; partial results may still be returned for other CIDs.
- Logs important events and warnings:
    - Missing molecules in SDF payloads
    - CID mismatches between JSON metadata and SDF properties
    - Conformer IDs being processed

##### Notes, constraints and assumptions

- This function uses an internal URL template to request PubChem similarity
    - results, then parses the response according to `output`:
    - SDF: parsed via SDFLoader.readDataFromText
    - JSON: parsed via json.loads and expected keys vary with `operation`
    - TXT: expected to contain one CID per line

- After extracting CIDs, the function calls get_compound(...) to obtain
    detailed compound (and optional conformer) data.
- Progress reporting, error continuation, and exact returned ALNPCompound
    structure are delegated to the underlying fetch/get routines.

- others
  - input 可选值为cid,smiles,InChI.

  - ouput 可选值为SDF,JSON或TXT

  - operation 可选值为record,cids,sids

  - others
| Option     | Type    | Meaning                                             | Default   |
| ---------- | ------- | --------------------------------------------------- | --------- |
| Threshold  | integer | minimum Tanimoto score for a hit                    | 90        |
| MaxSeconds | integer | maximum search time in seconds                      | unlimited |
| MaxRecords | integer | maximum number of hits                              | 2M        |
| listkey    | string  | restrict to matches within hits from a prior search | none      |

##### Example

- Simple usage:

```python

    cid = 2244
    compounds = get_similarity_compound(
                                input=EInputType.CID, 
                                value=cid, 
                                operation=EOperationType.CIDS, 
                                output=EOutputType.TXT)
    # compounds is a dict mapping "2244" -> ALNPCompound(...)

```

or

```python

    smiles = "CCCCCC1C(C(OC(=O)C(C(OC1=O)C)NC(=O)C2=C(C(=CC=C2)NC=O)O)C)OC(=O)CC(C)C"
    print(get_similarity_compound(input=EInputType.SMILES, 
                                  value=smiles, 
                                  operation=EOperationType.CIDS, 
                                  output=EOutputType.TXT))
    # compounds is a dict mapping smiles -> ALNPCompound(...)

```

- Fetch compounds with conformers:

```python

    smiles = "CCCCCC1C(C(OC(=O)C(C(OC1=O)C)NC(=O)C2=C(C(=CC=C2)NC=O)O)C)OC(=O)CC(C)C"
    print(get_similarity_compound(input=EInputType.CID, 
                                  value=cid_list[-1], 
                                  operation=EOperationType.CIDS, 
                                  output=EOutputType.TXT))
    # compounds is a dict mapping smiles -> ALNPCompound(...)

```

##### Types referenced

- ALNPCompound: container/dataclass representing a PubChem compound record and any attached conformer info.
- ALNPConformer: container/dataclass representing a single conformer (ID, parent CID, raw SDF).
Security and privacy
- Requests are made to the public PubChem REST endpoints; no credentials are required.
- Raw SDF content is stored in returned objects' ROW fields and may contain structural or identifier data;
    treat returned data according to your privacy/security policies.
>>>>>>> 307a51c7d6a82f56d6b8dda8235a1c23243dfa22
