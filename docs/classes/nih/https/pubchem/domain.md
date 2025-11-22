# PubChem Domain Types

`pybiotech.classes.nih.https.pubchem.domain` 定义了 PubChem PUG REST 输入路径中常用的 domain/namespace 子结构，包括结构搜索、快速搜索、xref、mass、fragment 等的枚举与校验型 dataclass，方便将 URL 片段映射为对象并进行反向拼接。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `StructureSearchType` / `StructureSearchSubType` | Enum | 控制结构搜索一级/二级关键词，如 `substructure/inchi`。 |
| `StructureSearch` | `BaseModel` | 将 `prefix/subtype` 片段映射为 `StructureSearchType` + `StructureSearchSubType`，并提供 `from_string`/`__str__`。 |
| `FastSearchType` / `FastSearchSubType` | Enum | 对应 `fastsuperstructure/smiles` 等快速检索前缀与子类型。 |
| `FastSearch` | `BaseModel` | 支持 `fastformula` 和 `prefix/subtype` 组合，保证语法合法。 |
| `XrefSubType` / `Xref` | Enum + `BaseModel` | `xref/<type>` 片段的解析器（例如 `xref/PubMedID`）。 |
| `MassType` / `MassMode` | Enum | 支持 `molecular_weight`、`equals/range` 等两个层级。 |
| `Mass` | `BaseModel` | 表示 `mass/type/mode/value1(/value2)`，支持 `range` 时 `value_2` 必填。 |
| `FASTFORMULA` | `str` | 约定常量 `"fastformula"`，简化特殊 fastsearch 解析。 |

---

### `StructureSearch`

**功能说明**

- 解析 `substructure/inchi` 形式的结构搜索参数，将前缀/子类型映射到预设的枚举。

**参数**

- `prefix`：`StructureSearchType`，如 `identity`/`similarity`。
- `subtype`：`StructureSearchSubType`，如 `smiles`/`inchi`。

**返回类型**

- `StructureSearch` 实例，`__str__` 会复原为 `prefix/subtype` 的字符串格式。

**异常**

- `ValueError`：`from_string` 输入不符合 `<prefix>/<subtype>` 或值不在枚举内时抛出。

**注意事项**

- ✅ 枚举确保拼接的子类型在 PubChem 支持范围内。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.domain import StructureSearch

ss = StructureSearch.from_string("superstructure/inchi")
print(str(ss))  # superstructure/inchi
```

---

### `FastSearch` & `FASTFORMULA`

**功能说明**

- `FastSearch` 表示 PubChem 快速结构检索片段，包括 `fastidentity/smiles` 以及特殊常量 `fastformula`（由 `FASTFORMULA` 提供）。

**参数**

- `prefix`：`FastSearchType` 枚举（可选）。
- `subtype`：`FastSearchSubType` 枚举（可选）。
- `fastformula`：布尔，`True` 时表示跳过 prefix/subtype。

**返回类型**

- `FastSearch` 实例；`__str__` 会按 `fastformula` 或 `prefix/subtype` 自动输出。

**异常**

- `ValueError`：`from_string` 无法识别传入的片段时抛出。

**注意事项**

- ⚠️ `fastformula` 与 `prefix/subtype` 互斥；模块内部在 `from_string("fastformula")` 直接返回特殊标志。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.domain import FastSearch

print(FastSearch.from_string("fastsuperstructure/smiles"))
print(FastSearch.from_string("fastformula"))
```

---

### `Xref`

**功能说明**

- 解析 `xref/<tag>` 片段，确保只允许官方列出的 `XrefSubType`，并提供 `__str__` 供 URL 拼接。

**参数**

- `subtype`：`XrefSubType` 枚举（例如 `PubMedID`、`GeneID`）。

**返回类型**

- `Xref` 实例，`__str__` 输出 `xref/<subtype>`。

**异常**

- `ValueError`：字符串格式错误或 subtype 不合法。

**注意事项**

- ✅ 所有 subtype 通过枚举校验，避免拼写错误。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.domain import Xref

print(Xref.from_string("xref/PubMedID"))
```

---

### `Mass`

**功能说明**

- 表示 `mass/<type>/<mode>/value` 或 `mass/<type>/range/value1/value2` 的片段，在 `range` 时强制要求 `value_2`。

**参数**

- `mass_type`：`MassType`，如 `molecular_weight`。
- `mode`：`MassMode`（`equals`/`range`）。
- `value_1`：主数值（`float`）。
- `value_2`：`range` 情况下的上限，`equals` 时必须为 `None`。

**返回类型**

- `Mass` 实例；`__str__` 输出完整路径。

**异常**

- `ValueError`：`value_2` 与 `mode` 不匹配或 `from_string` 解析失败。

**注意事项**

- ⚠️ `field_validator` 保证 `range` 有上界、`equals` 无上界。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.domain import Mass

print(Mass.from_string("exact_mass/range/100/200"))
```
